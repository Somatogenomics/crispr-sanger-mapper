library(shiny)

for (f in c("R/utils.R","R/abif.R","R/qc.R","R/alignment.R","R/analysis.R","R/plots.R")) source(f, local=FALSE)

ui <- fluidPage(
  tags$head(tags$link(rel="stylesheet", type="text/css", href="styles.css")),
  titlePanel("CRISPR Sanger Mapper"),
  div(class="subtitle", "Browser-local AB1 QC, reference mapping, sgRNA/PAM validation, and cut-site trace screening"),
  sidebarLayout(
    sidebarPanel(width=3,
      fileInput("ref", "Reference FASTA", accept=c(".fa",".fasta",".fna","text/plain")),
      uiOutput("ref_record_ui"),
      fileInput("ab1", "Sanger chromatograms (.ab1)", accept=c(".ab1"), multiple=TRUE),
      textAreaInput("guides", "sgRNA protospacer(s), one per line", rows=4,
                    placeholder="ACGTTGCTACGATCCTGATG"),
      textInput("pam", "PAM motif (IUPAC)", value="NGG"),
      tags$hr(),
      numericInput("mott", "Mott error-probability cutoff", value=0.05, min=.001, max=.2, step=.01),
      numericInput("minq", "Fallback minimum Phred", value=20, min=0, max=60),
      numericInput("minlen", "Minimum HQ window (bp)", value=20, min=10, max=200),
      numericInput("het", "Secondary / primary peak threshold", value=.30, min=.05, max=.95, step=.05),
      actionButton("run", "Analyze", class="btn-primary btn-block"),
      tags$hr(),
      p(class="privacy", strong("Privacy:"), "uploaded files are processed in your browser when deployed with Shinylive; they are not sent to an analysis server."),
      p(class="note", "Automated calls are screening/QC outputs. Poor or phase-shifted Sanger traces require manual review or amplicon deep sequencing.")
    ),
    mainPanel(width=9,
      tabsetPanel(
        tabPanel("Batch summary",
          br(), uiOutput("run_note"),
          div(class="table-wrap", tableOutput("summary")),
          downloadButton("download_summary","Download summary CSV")
        ),
        tabPanel("Sample review",
          br(), fluidRow(column(5,selectInput("selected","Sample × guide", choices=NULL)),
                         column(7,uiOutput("call_box"))),
          uiOutput("site_evidence"),
          plotOutput("trace", height="430px"),
          h4("Alignment"), pre(class="alignment", textOutput("alignment")),
          h4("Reference mismatches"), tableOutput("mismatches")
        ),
        tabPanel("Method",
          br(),
          h4("Analysis stages"),
          tags$ol(
            tags$li("Parse ABIF tags: PBAS, PCON, PLOC, FWO_, and analyzed trace channels DATA9–DATA12."),
            tags$li("Quality trim using Mott scoring, with a longest-Q fallback."),
            tags$li("Align both read orientations to the selected FASTA record."),
            tags$li("Require a strand-correct PAM adjacent to the supplied guide before deriving a SpCas9 cut coordinate."),
            tags$li("Map BOTH reference bases flanking the predicted Cas9 cleavage boundary into the original chromatogram, accounting for read and sgRNA strand."),
            tags$li("Assess the pre-cut baseline (-35 to -9 bases), onset near the cut (-5 to +8), and sustained post-cut signal (+3 to +28), all relative to the sequencing direction and restricted to mapped bases."),
            tags$li("Require a clean-enough pre-cut baseline and a significant increase in sustained secondary/primary ratio after the validated cut. Diffuse high peaks elsewhere never trigger the cut-localized call."),
            tags$li("Separate low-quality or unbracketed-cut results from possible editing; report nearby alignment gaps independently for manual review.")
          ),
          h4("Interpretation"),
          p("WT-like means no automated indel signature was detected in a read that passed QC. It does not prove the absence of a low-frequency edit."),
          p("A cut-localized phase change is a screening clue, not a diagnosis of heterozygosity or mosaicism. Without a matched WT chromatogram, the app does not subtract background or estimate allele frequencies."),
          p("Reference mismatches away from the cut should be checked in unedited strain-matched controls before being attributed to editing.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  refs <- reactive({
    req(input$ref)
    read_fasta_base(input$ref$datapath)
  })
  output$ref_record_ui <- renderUI({
    r <- refs(); choices <- setNames(seq_along(r), vapply(r,function(z) paste0(z$name," [",nchar(z$seq)," bp]"),character(1)))
    selectInput("ref_record", "FASTA record", choices=choices, selected=1)
  })

  analysis <- eventReactive(input$run, {
    req(input$ref, input$ab1, input$ref_record)
    gs <- parse_guides(input$guides)
    validate(need(length(gs)>0,"Enter at least one sgRNA protospacer."))
    r <- refs()[[as.integer(input$ref_record)]]
    withProgress(message="Analyzing chromatograms", value=0, {
      ans <- analyse_batch(input$ab1$datapath, r$name, r$seq, gs, toupper(input$pam),
                           sample_names=input$ab1$name,
                           mott_cutoff=input$mott, min_hq_len=input$minlen,
                           min_qual=input$minq, het_threshold=input$het)
      incProgress(1)
      ans
    })
  })

  observeEvent(analysis(), {
    a <- analysis(); ch <- names(a$results)
    updateSelectInput(session,"selected",choices=setNames(ch, gsub("\\|\\|"," — ",ch)), selected=ch[1] %||% character())
  })

  output$run_note <- renderUI({
    a <- analysis()
    if (!length(a$errors)) return(div(class="ok-box", "Analysis completed without parser/runtime errors."))
    div(class="warn-box", strong("Some files could not be analyzed:"), tags$ul(lapply(names(a$errors),function(n) tags$li(n,": ",a$errors[[n]]))))
  })

  output$summary <- renderTable({ analysis()$summary }, striped=TRUE, bordered=TRUE, spacing="xs", digits=3)
  output$download_summary <- downloadHandler(
    filename=function() paste0("crispr_sanger_mapper_summary_",Sys.Date(),".csv"),
    content=function(file) write.csv(analysis()$summary,file,row.names=FALSE)
  )

  selected_result <- reactive({ req(input$selected); analysis()$results[[input$selected]] })
  output$call_box <- renderUI({
    x <- selected_result()$summary[1,]
    cls <- if (x$status=="FAIL") "warn-box" else if (grepl("POSSIBLE",x$call)) "candidate-box" else "ok-box"
    div(class=cls,
      strong(x$call), tags$br(),
      sprintf("Status: %s | strand: %s | identity: %s | mean Q: %s | cut covered: %s",
              x$status, x$read_strand, ifelse(is.na(x$identity),"NA",sprintf("%.1f%%",100*x$identity)), x$mean_q, x$cut_covered),
      if (nzchar(x$manual_review)) tagList(tags$br(),em(x$manual_review))
    )
  })
  output$site_evidence <- renderUI({
    a <- selected_result(); x <- a$summary[1,]
    if (is.null(a$guide_hit)) return(div(class="warn-box", "No uniquely PAM-validated sgRNA was found; cut-centered inference is disabled."))
    e <- a$cut_evidence
    fmt <- function(z, digits=2) if (length(z)==0 || !is.finite(z)) "NA" else sprintf(paste0("%.",digits,"f"),z)
    div(class="note", tags$strong("Guide-anchored evidence"), tags$br(),
      sprintf("sgRNA reference %d–%d (%s strand); PAM %s at %d–%d; predicted SpCas9 cleavage between reference %d|%d.",
              a$guide_hit$guide_start,a$guide_hit$guide_end,a$guide_hit$strand,
              a$guide_hit$pam,a$guide_hit$pam_start,a$guide_hit$pam_end,
              a$guide_hit$cut_after,a$guide_hit$cut_after+1L),
      tags$br(),
      if (!isTRUE(x$cut_covered)) "Both cut-flanking bases are not adjacent in the mapped read; no mixed-indel inference is permitted." else
        sprintf("Aligned cut baseline: n=%d median ratio=%s; near-cut mixed bases (-5 to +8): %d; downstream (+3 to +28): %d/%d; median ratio=%s; increase=%s. Cut-localized phase pattern: %s.",
                e$baseline_n,fmt(e$baseline_median),e$onset_n,e$post_double,e$post_n,
                fmt(e$post_median),fmt(e$delta),if(e$pattern) "YES (candidate)" else "NO"),
      tags$br(), tags$small("All offsets follow original sequencing direction. Secondary/primary ratios are descriptive signals; results are not confirmed genotypes.")
    )
  })
  output$trace <- renderPlot({ plot_trace(selected_result()) }, res=110)
  output$alignment <- renderText({ alignment_strings(selected_result()) })
  output$mismatches <- renderTable({
    m <- selected_result()$mismatches
    if (!nrow(m)) data.frame(message="No aligned base mismatches detected") else m
  }, striped=TRUE, bordered=TRUE, spacing="xs")
}

shinyApp(ui, server)
