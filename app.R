library(shiny)

# Chromium + Shinylive workaround: standard downloadButton's HTML `download`
# attribute may break file downloads in Chrome/Edge when running via webR.
# See https://forum.posit.co/t/shinylive-download-works-with-firefox-fails-on-chrome/198542
downloadButton <- function(...) {
  button <- shiny::downloadButton(...)
  button$attribs$download <- NULL
  button
}

for (f in c("R/utils.R","R/abif.R","R/qc.R","R/alignment.R","R/analysis.R","R/plots.R","R/reporting.R")) source(f, local=FALSE)

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
          h4("Screening and QC summary"),
          div(class="table-wrap", tableOutput("summary_short")),
          div(class="download-actions",
              downloadButton("download_summary", "Complete results CSV"),
              downloadButton("download_summary_tsv", "Complete results TSV"),
              downloadButton("download_parameters", "Analysis settings CSV"),
              downloadButton("download_errors", "Parser errors CSV")),
          h4("Batch quality figure"),
          plotOutput("batch_quality", height="390px"),
          downloadButton("download_batch_png", "Download batch QC PNG"),
          tags$details(tags$summary("Show all result columns"),
                       div(class="table-wrap", tableOutput("summary")))
        ),
        tabPanel("Sample review",
          br(), fluidRow(column(5,selectInput("selected","Sample × guide", choices=NULL)),
                         column(7,uiOutput("call_box"))),
          div(class="download-actions",
              downloadButton("download_sample_pdf", "Complete sample report PDF"),
              downloadButton("download_sample_csv", "Sample evidence CSV"),
              downloadButton("download_per_base", "Per-base metrics TSV")),
          h4("1. Chromatogram around the PAM-validated cut"),
          plotOutput("trace", height="420px"),
          downloadButton("download_trace_png", "Download chromatogram PNG"),
          h4("2. Peak mixing: onset and sustained signal relative to cut"),
          plotOutput("peak_ratio", height="380px"),
          downloadButton("download_ratio_png", "Download peak-ratio PNG"),
          h4("3. Phred quality, HQ window, and alignment coverage"),
          plotOutput("quality_profile", height="360px"),
          downloadButton("download_quality_png", "Download quality PNG"),
          h4("Alignment"), pre(class="alignment", textOutput("alignment")),
          h4("Cut-site evidence"), div(class="table-wrap",tableOutput("evidence")),
          h4("Alignment gaps near cut"), tableOutput("gaps"),
          h4("Reference mismatches"), tableOutput("mismatches"),
          p(class="note", "Figures display the original sequencing direction. An approximate/unresolved cut is marked for visual review only. QC failures must not be interpreted as a negative genotype.")
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
          h4("Downloadable outputs"),
          p("The batch summary exports CSV/TSV plus analysis settings; each sample exports an evidence CSV, base-by-base TSV with alignment and peak metrics, three PNG figures, and a multi-page PDF report. Downloads are produced in browser-local webR memory when deployed as Shinylive."),
          p("Chrome/Edge compatibility: download buttons omit the HTML download attribute to avoid a known Shinylive browser bug."),
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
      # Freeze the settings used at Analyze time. Plot overlays, PDF reports,
      # and downloads must describe the SAME thresholds as the genotype screen.
      ans$parameters <- list(mott_cutoff=input$mott,
                             fallback_min_phred=input$minq,
                             min_hq_window_bp=input$minlen,
                             secondary_primary_threshold=input$het,
                             min_mean_hq_phred=input$minmeanq,
                             aln_match=input$aln_match,
                             aln_mismatch=input$aln_mismatch,
                             aln_gap_open=input$aln_open,
                             aln_gap_extend=input$aln_extend)
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

  output$summary_short <- renderTable({
    x <- analysis()$summary
    if (!nrow(x)) return(data.frame(message="No completed analyses; inspect parser errors."))
    x[, c("sample","guide","status","call","strand","hq_len","mean_q",
          "identity","cut_mapping","pre_scored_n","post_scored_n",
          "cut_localized_pattern"), drop=FALSE]
  }, striped=TRUE, bordered=TRUE, spacing="xs", digits=3)
  output$summary <- renderTable({ analysis()$summary }, striped=TRUE, bordered=TRUE, spacing="xs", digits=3)
  output$download_summary <- downloadHandler(
    filename=function() paste0("crispr_sanger_mapper_summary_",Sys.Date(),".csv"),
    content=function(file) utils::write.csv(analysis()$summary,file,row.names=FALSE,na=""),
    contentType="text/csv"
  )
  output$download_summary_tsv <- downloadHandler(
    filename=function() paste0("crispr_sanger_mapper_summary_",Sys.Date(),".tsv"),
    content=function(file) utils::write.table(analysis()$summary,file,sep="\t",row.names=FALSE,quote=TRUE,na=""),
    contentType="text/tab-separated-values"
  )
  output$download_parameters <- downloadHandler(
    filename=function() paste0("crispr_sanger_mapper_settings_",Sys.Date(),".csv"),
    content=function(file) write_settings_csv(analysis()$parameters,file),
    contentType="text/csv"
  )
  output$download_errors <- downloadHandler(
    filename=function() paste0("crispr_sanger_mapper_parser_errors_",Sys.Date(),".csv"),
    content=function(file) {
      err <- analysis()$errors
      d <- data.frame(sample=names(err),error=unlist(err, use.names=FALSE),stringsAsFactors=FALSE)
      utils::write.csv(d,file,row.names=FALSE)
    }, contentType="text/csv"
  )
  output$batch_quality <- renderPlot({ plot_batch_quality(analysis()$summary) }, res=115)
  output$download_batch_png <- downloadHandler(
    filename=function() paste0("crispr_sanger_batch_QC_",Sys.Date(),".png"),
    content=function(file) save_figure_png(file,function() plot_batch_quality(analysis()$summary)),
    contentType="image/png"
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
  # Never silently update figure thresholds before the analysis is re-run.
  plotted_threshold <- reactive(analysis()$parameters$secondary_primary_threshold)
  output$trace <- renderPlot({ plot_trace(selected_result(), het_threshold=plotted_threshold()) }, res=110)
  output$peak_ratio <- renderPlot({ plot_peak_ratio(selected_result(), het_threshold=plotted_threshold()) }, res=110)
  output$quality_profile <- renderPlot({ plot_quality_profile(selected_result()) }, res=110)
  output$download_trace_png <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"png","trace"),
    content=function(file) save_figure_png(file,function() plot_trace(selected_result(),het_threshold=plotted_threshold())),
    contentType="image/png"
  )
  output$download_ratio_png <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"png","cut_peak_ratio"),
    content=function(file) save_figure_png(file,function() plot_peak_ratio(selected_result(),het_threshold=plotted_threshold())),
    contentType="image/png"
  )
  output$download_quality_png <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"png","phred_quality"),
    content=function(file) save_figure_png(file,function() plot_quality_profile(selected_result())),
    contentType="image/png"
  )
  output$download_sample_pdf <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"pdf","screening_report"),
    content=function(file) write_sample_report_pdf(selected_result(),file,het_threshold=plotted_threshold(),
                                                    settings=analysis()$parameters),
    contentType="application/pdf"
  )
  output$download_sample_csv <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"csv","evidence"),
    content=function(file) write_sample_evidence(selected_result(),file),
    contentType="text/csv"
  )
  output$download_per_base <- downloadHandler(
    filename=function() sample_export_name(selected_result(),"tsv","per_base"),
    content=function(file) write_sample_per_base(selected_result(),file),
    contentType="text/tab-separated-values"
  )
  # Keep downloads active after initial render even if another tab is selected.
  for (id in c("download_summary","download_summary_tsv","download_parameters","download_errors",
               "download_batch_png","download_trace_png","download_ratio_png",
               "download_quality_png","download_sample_pdf","download_sample_csv",
               "download_per_base")) {
    outputOptions(output, id, suspendWhenHidden=FALSE)
  }
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
