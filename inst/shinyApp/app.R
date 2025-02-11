# Load required libraries
suppressPackageStartupMessages({
  
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(CMplot)
  library(DT)
  library(ggrepel)
  
})

rm(list = ls())
gc()

# Load functions
package_extdata_path <- paste(system.file(package = "BRAINscape"), "extdata", sep = "/")

ui <- fluidPage(
  titlePanel("Analysis Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ethnicity", "Select Ethnicity:", choices = c("Hispanics", "NonHispanics")),
      selectInput("analysisType", "Select Analysis Type:", 
                  choices = c("eQTL", "Differential Expression Analysis")),
      
      # Conditional UI for eQTL analysis
      conditionalPanel(
        condition = "input.analysisType == 'eQTL'",
        radioButtons("eqtlType", "Select eQTL Type:", choices = c("cis", "trans")),
        selectizeInput("gene", "Select Gene (eQTL):", choices = NULL, multiple = FALSE),
        sliderInput("pval_nominal", "Nominal P-Value Threshold:", min = 0, max = 1, value = 0.05, step = 0.01),
        actionButton("filterBtn", "Apply Filter")
      ),
      
      # Conditional UI for DESeq2 analysis
      conditionalPanel(
        condition = "input.analysisType == 'Differential Expression Analysis'",
        selectInput("p_adjusted", label = "Cutoff for p-adjusted:", choices = c("0.00001" = 0.00001,"0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0)),
        numericInput("base_mean", label = "Minimal base mean:", value = 0),
        numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary Table", DTOutput("summaryTable")),
        tabPanel("Manhattan Plot", plotOutput("manhattanPlot")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive variables for data loading
  eqtl_data <- reactive({
    req(input$analysisType == "eQTL")
    load_eqtl(package_extdata_path, input$ethnicity, input$eqtlType)
  })
  
  deseq2_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
    withProgress(message = "Loading DESeq2 data...", value = 0, {
      load_deseq2(package_extdata_path, input$ethnicity)
    })
  })
  
  expression_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
    withProgress(message = "Loading Expression data...", value = 0, {
      load_expression_data(package_extdata_path, input$ethnicity)
    })
  })
  
  pheno_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
      load_pheno_data(package_extdata_path)
  })
  
  sample_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
    load_samples(package_extdata_path, input$ethnicity)
  })
  
  # Observing and updating gene dropdown for eQTL analysis
  observe({
    if (input$analysisType == "eQTL") {
      req(eqtl_data(), input$pval_nominal)
      filtered_genes <- eqtl_data() %>%
        filter(pval < input$pval_nominal) %>%
        pull(gene) %>%
        unique()
      updateSelectizeInput(session, "gene", choices = filtered_genes, server = TRUE)
    } else {
      updateSelectizeInput(session, "gene", choices = NULL, server = TRUE)
    }
  })
  
  # Filtering data based on analysis type
  filtered_data <- reactive({
    if (input$analysisType == "eQTL") {
      req(input$gene, input$pval_nominal)
      eqtl_data() %>% filter(gene == input$gene & pval <= input$pval_nominal)
    } else {
      req(input$p_adjusted, input$base_mean, input$log2fc)
      deseq2_data() %>% filter(padj <= as.numeric(input$p_adjusted) & baseMean >= input$base_mean & abs(log2FoldChange) >= input$log2fc)
    }
  })
  
  # Summary table
  output$summaryTable <- renderDT({
    req(filtered_data())
    datatable(filtered_data())
  })
  
  # Volcano Plot
  output$volcanoPlot <- renderPlot({
    req(input$analysisType == "Differential Expression Analysis", filtered_data())
    res <- filtered_data()
    res$diffexpressed <- ifelse(res$log2FoldChange > 0, "UP", "DOWN")
    colnames(res)[1] <- "gene_symbol"
    res$delabel <- ifelse(res$gene_symbol %in% head(res[order(res$padj), "gene_symbol"], 30), res$gene_symbol, NA)
    ggplot(data = res, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
      geom_point(size = 2) +
      # Updated colors for UP, DOWN, and NO genes
      scale_color_manual(values = c("DOWN" = "#00AFBB", "UP" = "#bb0c00"), 
                         labels = c("Downregulated",  "Upregulated")) + 
      coord_cartesian(ylim = c(0, 15), xlim = c(-2, 2)) + 
      labs(color = 'Expression Status', # Updated legend title
           x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
      scale_x_continuous(breaks = seq(-2, 2, 0.5)) + 
      ggtitle(paste0('Volcano Plot for ', input$ethnicity, ' ethnicity')) +
      geom_text_repel(max.overlaps = 20, size = 5, box.padding = 0.5, max.time = 3) +  # Improve readability +
      theme_minimal(base_size = 14)  # Better font clarity
    }, height = 700, width = 700)
  
  
  # Manhattan Plot
  output$manhattanPlot <- renderPlot({
    req(filtered_data())
    
    if (input$analysisType == "eQTL") {
      
      # Manhattan plot for eQTL
      manhattan_df <- data.frame(
        SNP = filtered_data()$snp,
        Chromosome = filtered_data()$CHR,
        Position = filtered_data()$POS,
        P = filtered_data()$pval
      )
      
      if (nrow(manhattan_df) <= 1) {
        output$manhattanPlotMessage <- renderText("The Manhattan plot can't be generated for less than one SNP.")
        plot(1, type = "n", xlab = "", ylab = "", main = "The Manhattan plot can't be generated for less than one SNP.")
        return()
      }
      
      output$manhattanPlotMessage <- renderText("")  # Clear any previous message
      
      manhattan_df <- manhattan_df %>% arrange(P)
      
      snp_annot <- manhattan_df$SNP[1:5]
      CMplot(manhattan_df,
             type="p",
             plot.type="m",
             col = c("blue", "red"),
             LOG10=TRUE,
             file = NULL,
             cex = 0.6, main = paste0("Manhattan Plot of Cis-eQTL Analysis for ", input$gene, " Gene in ", input$ethnicity, " Cohort"),
             file.output=F,verbose=F,
             width=14,height=6,chr.labels.angle=45,
             highlight=snp_annot, 
             highlight.text=snp_annot, 
             highlight.text.cex=1.4,
             highlight.col = "red")
      
      
    }
  })

}


shinyApp(ui = ui, server = server)
