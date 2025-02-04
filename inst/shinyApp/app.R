# Load required libraries
suppressPackageStartupMessages({
  
  library(shiny)
  library(ggplot2)
  library(ComplexHeatmap)
  library(DT)
  library(DESeq2)
  library(circlize)
  library(dplyr)
  library(pheatmap)  # For heatmap visualization
  library(CMplot)
  
})

rm(list = ls())
gc()

# Load functions
# source(system.file("R/load_data.R", package = "BRAINscape"))
source("/mnt/vast/hpc/homes/vr2592/BRAIN/data_analysis/BRAINscape/R/load_data.R")
package_extdata_path <- paste("/mnt/vast/hpc/homes/vr2592/BRAIN/data_analysis/BRAINscape", "inst", "extdata", sep = "/")

ui <- fluidPage(
  titlePanel("Analysis Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ethnicity", "Select Ethnicity:", choices = c("Hispanics", "NonHispanics")),
      selectInput("analysisType", "Select Analysis Type:", 
                  choices = c("eQTL", "GWAS", "Differential Expression Analysis")),
      
      # Conditional UI for eQTL analysis
      conditionalPanel(
        condition = "input.analysisType == 'eQTL'",
        radioButtons("eqtlType", "Select eQTL Type:", choices = c("cis", "trans")),
        selectInput("gene", "Select Gene (eQTL):", choices = NULL),
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
        tabPanel("QQ Plot", plotOutput("qqPlot")),
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
  
  gwas_data <- reactive({
    req(input$analysisType == "GWAS")
    withProgress(message = "Loading GWAS data...", value = 0, {
      # Step 1: Initialize
      incProgress(0.2, detail = "Initializing data loading...")
      Sys.sleep(0.25) # Simulate some processing time
      
      # Step 2: Load the data
      gwas <- load_gwas(package_extdata_path, input$ethnicity)
      incProgress(0.5, detail = "Loading GWAS data into memory...")
      Sys.sleep(0.25) # Simulate some processing time
      
      # Step 3: Finalizing
      incProgress(0.3, detail = "Finalizing...")
      Sys.sleep(0.25) # Simulate some processing time
      
      # Return the loaded data
      gwas
    })
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
      updateSelectInput(session, "gene", choices = filtered_genes)
    } else {
      updateSelectInput(session, "gene", choices = NULL)
    }
  })
  
  # Filtering data based on analysis type
  filtered_data <- reactive({
    if (input$analysisType == "eQTL") {
      req(input$gene, input$pval_nominal)
      eqtl_data() %>% filter(gene == input$gene & pval <= input$pval_nominal)
    } else if (input$analysisType == "GWAS") {
      gwas_data()
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
      coord_cartesian(ylim = c(0, 15), xlim = c(-8, 8)) + 
      labs(color = 'Expression Status', # Updated legend title
           x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
      scale_x_continuous(breaks = seq(-8, 8, 2)) + 
      ggtitle(paste0('Volcano Plot for ', input$ethnicity, ' ethnicity')) +
      geom_text_repel(max.overlaps = Inf, size = 5)
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
      
      
    } else if (input$analysisType == "GWAS") {
      
      # Manhattan plot for GWAS
      req(input$analysisType == "GWAS", filtered_data())
      
      # Prepare GWAS data
      gwas_data <- filtered_data()
      manhattan_df <- data.frame(
        SNP = gwas_data$snp,
        Chromosome = as.factor(gwas_data$chr),
        Position = gwas_data$pos,
        P = gwas_data$`p-value`
      )
      
      manhattan_df <- manhattan_df %>% arrange(P)
      
      message("Manhattan data is prepared and the dimension is = ", dim(manhattan_df))
      
      manhattan_df <- head(manhattan_df, n = 50000)
      
      message("Now the dimension is = ", dim(manhattan_df))
      
      if (nrow(manhattan_df) == 0) {
        plot(1, type = "n", main = "No data available for the selected criteria.", xlab = "", ylab = "")
        return()
      }
      
      # Calculate significance thresholds
      significant_threshold <- 0.05 / nrow(manhattan_df)
      suggestive_threshold <- significant_threshold * 100
      SNPs <- manhattan_df$SNP[1:5]
      
      manhattan_df <- manhattan_df %>%
        mutate(Chromosome = as.numeric(gsub("chr", "", Chromosome))) %>%
        arrange(Chromosome, Position)
      
      CMplot(manhattan_df,
             type="p",
             plot.type="m",
             col = c("blue", "red"),
             LOG10=TRUE,
             file = NULL,
             cex = 0.6,
             main = paste0("Genome-wide Association Study (GWAS) Manhattan Plot of top 50K hits for ", input$ethnicity, " Cohort"),
             file.output=F,
             verbose=F,
             width=14,
             height=6,
             chr.labels.angle=45,
             highlight=SNPs,
             highlight.text=SNPs,
             highlight.text.cex=1.4,
             highlight.col = "red")
      
      
      # 
      # # Combine plots using gridExtra
      # gridExtra::grid.arrange(manhattan_plot, qq_plot(), nrow = 2)
    }
  })
  
  output$qqPlot <- renderPlot({
    req(filtered_data())
    
    if (input$analysisType == "GWAS") {
      # Step 1: Start processing
      output$qqPlotMessage <- renderUI({
        HTML("<b>Step 1:</b> Estimating chi distribution for the provided data...")
      })
      Sys.sleep(0.5) # Simulate some processing delay
      
      # Step 2: Estimate chi distribution
      chi <- qchisq(1 - filtered_data()$`p-value`, 1)
      df <- data.frame(chi)
      
      # Step 3: Calculate lambda
      output$qqPlotMessage <- renderUI({
        HTML("<b>Step 2:</b> Calculating lambda value...")
      })
      Sys.sleep(0.5) # Simulate some processing delay
      lambda <- median(df$chi, na.rm = TRUE) / qchisq(0.5, 1)
      
      # Step 4: Prepare data for QQ plot
      output$qqPlotMessage <- renderUI({
        HTML("<b>Step 3:</b> Preparing data for QQ plot...")
      })
      Sys.sleep(0.5) # Simulate some processing delay
      observed <- sort(filtered_data()$`p-value`)
      lobs <- -log10(observed)
      
      expected <- -log10(ppoints(length(observed)))
      
      # Step 5: Render QQ plot
      output$qqPlotMessage <- renderUI({
        HTML("<b>Step 4:</b> Generating QQ plot...")
      })
      Sys.sleep(1) # Simulate some processing delay
      
      plot(
        c(0, max(expected)), c(0, max(lobs)), col = "red", lwd = 3, type = "l",
        xlab = "Expected (-logP)", ylab = "Observed (-logP)", 
        xlim = c(0, max(expected)), ylim = c(0, max(lobs)), las = 1, xaxs = "i", yaxs = "i", bty = "l"
      )
      points(expected, lobs, pch = 23, cex = 0.4, bg = "black")
      text(1, max(lobs) * 0.9, paste0("Lambda = ", round(lambda, 3)), adj = 0, cex = 0.8, font = 2, col = "blue")
      
      # Step 6: Completion
      output$qqPlotMessage <- renderUI({
        HTML("<b>Step 5:</b> QQ plot successfully generated.")
      })
    }
  })
  
  # UI to display messages
  output$qqPlotMessage <- renderUI({
    HTML("<b>Initializing QQ plot generation...</b>")
  })
  
}


shinyApp(ui = ui, server = server)
