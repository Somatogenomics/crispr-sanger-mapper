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
                    placeholder="20-nt protospacer (without PAM)"),
      textInput("pam", "PAM motif (IUPAC)", value="NGG"),
      tags$hr(),
      numericInput("mott", "Mott error-probability cutoff", value=0.05, min=.001, max=.2, step=.01),
      numericInput("minq", "Fallback minimum Phred", value=20, min=0, max=60),
      numericInput("minlen", "Minimum HQ window (bp)", value=20, min=10, max=200),
      numericInput("het", "Secondary/primary peak-height ratio", value=.30, min=.05, max=.95, step=.05),
      numericInput("minmeanq", "Minimum mean quality (HQ window)", value=15, min=0, max=60, step=1),
      tags$details(
        tags$summary("Alignment scoring"),
        numericInput("aln_match", "Match", value=2, step=1),
        numericInput("aln_mismatch", "Mismatch", value=-2, step=1),
        numericInput("aln_open", "Gap open (first gap position)", value=-5, step=1),
        numericInput("aln_extend", "Gap extend (each further position)", value=-1, step=1)
      ),
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
          plotOutput("trace", height="390px"),
          h4("Alignment"), pre(class="alignment", textOutput("alignment")),
          h4("Cut-site evidence"), tableOutput("evidence"),
          h4("Alignment gaps near cut"), tableOutput("gaps"),
          h4("Reference mismatches"), tableOutput("mismatches")
        ),
        tabPanel("Method",
          br(),
          h4("Analysis stages"),
          tags$ol(
            tags$li("Parse ABIF tags: PBAS, PCON, PLOC, FWO_, and analyzed trace channels DATA9–DATA12."),
            tags$li("Quality trim using Mott scoring, with a longest-Q fallback."),
            tags$li("Align both read orientations to the selected FASTA record (local alignment, affine gap scores)."),
            tags$li("Require a strand-correct PAM adjacent to the supplied guide before deriving a SpCas9 cut coordinate."),
            tags$li("Map BOTH reference bases immediately flanking the SpCas9 cut into original sequencing-direction coordinates. A one-sided or extrapolated cut can guide visual review but cannot support an automated mixed-indel or WT-like call."),
            tags$li("For mixed-phase calls, require a quality-filtered cleaner baseline before the cut, onset within -5..+8 bp, and sustained elevated secondary/primary signal in +3..+28 bp after the cut (in read direction)."),
            tags$li("A high-quality clean alignment gap bracketing the expected cut can yield a separate provisional clean-indel flag, even if the deleted base prevents mapping the immediate cut flanks. Diffuse or late-onset mixed signal never yields a cut-associated mixed-indel call."),
            tags$li("Keep inadequate-quality, inadequate-cut-window and ambiguous PAM cases as screening failures. Report reference mismatches for comparison with WT controls; do not automatically treat SNPs as edits.")
          ),
          h4("Interpretation"),
          p("WT-like means no automated indel signature was detected in an adequately mapped read; it does NOT prove absence of editing, especially low-frequency alleles or SNVs."),
          p("Possible heterozygous/mosaic indel requires PAM-validated cut-localized phase-change evidence and sufficient mapped baseline/onset/post-cut bases. It is a hypothesis, not proof of genotype."),
          p("Alignment gaps near cut are screened separately: their locations and allele identity require manual review and independent confirmation."),
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
    validate(need(input$aln_open <= 0 && input$aln_extend <= 0 && input$aln_match > 0,
                  "Gap scores must be <= 0 and the match score > 0."))
    r <- refs()[[as.integer(input$ref_record)]]
    withProgress(message="Analyzing chromatograms", value=0, {
      ans <- analyse_batch(input$ab1$datapath, r$name, r$seq, gs, toupper(input$pam),
                           sample_names=input$ab1$name,
                           mott_cutoff=input$mott, min_hq_len=input$minlen,
                           min_qual=input$minq, min_mean_q=input$minmeanq, het_threshold=input$het,
                           match=input$aln_match, mismatch=input$aln_mismatch,
                           gap_open=input$aln_open, gap_extend=input$aln_extend)
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
    cls <- if (x$status=="PASS" && grepl("WT_LIKE",x$call)) "ok-box" else if (x$status=="FAIL") "warn-box" else "candidate-box"
    div(class=cls,
      strong(x$call), tags$br(),
      sprintf("Status: %s | strand: %s | identity: %s | mean Q: %s | cut mapping: %s",
              x$status, x$strand, ifelse(is.na(x$identity),"NA",sprintf("%.1f%%",100*x$identity)), x$mean_q, x$cut_mapping),
      tags$br(), sprintf("Guide strand: %s | validated PAM: %s | cut: %s | near-cut onset: %s",
                        x$guide_strand, x$pam, x$cut_between, x$onset_first_offset),
      if (nzchar(x$manual_review)) tagList(tags$br(),em(x$manual_review))
    )
  })
  output$trace <- renderPlot({ plot_trace(selected_result(), het_threshold=input$het) }, res=110)
  output$alignment <- renderText({ alignment_strings(selected_result()) })
  output$evidence <- renderTable({
    x <- selected_result()$summary
    x[,c("cut_mapping","pre_scored_n","pre_median_ratio","onset_mixed_n",
         "onset_first_offset","post_scored_n","post_double_n","post_double_fraction",
         "post_median_ratio","post_minus_pre_ratio","cut_localized_pattern"),drop=FALSE]
  }, striped=TRUE, bordered=TRUE, spacing="xs", digits=3)
  output$gaps <- renderTable({
    g <- selected_result()$gaps
    if (!nrow(g)) data.frame(message="No alignment gaps near the cut") else g
  }, striped=TRUE, bordered=TRUE, spacing="xs")
  output$mismatches <- renderTable({
    m <- selected_result()$mismatches
    if (!nrow(m)) data.frame(message="No aligned base mismatches detected") else m
  }, striped=TRUE, bordered=TRUE, spacing="xs")
}

shinyApp(ui, server)
