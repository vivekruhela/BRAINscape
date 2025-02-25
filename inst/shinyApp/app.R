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
      uiOutput("ethnicity_ui"),
      selectInput("analysisType", "Select Analysis Type:", 
                  choices = c("eQTL", "Differential Expression Analysis")),
      uiOutput("ethnicity_ui"),
      
      # Conditional UI for eQTL analysis
      conditionalPanel(
        condition = "input.analysisType == 'eQTL'",
        radioButtons("eqtlType", "Select eQTL Type:", choices = c("cis", "trans")),
        selectizeInput("gene", "Select Gene (eQTL):", choices = NULL, multiple = FALSE),
        uiOutput("eqtlFilters")
      ),
      
      # Conditional UI for DESeq2 analysis
      conditionalPanel(
        condition = "input.analysisType == 'Differential Expression Analysis'",
        selectInput("p_adjusted", label = "Adjusted P-Value Threshold:", 
                    choices = c("0.00001" = 0.00001,"0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0)),
        numericInput("base_mean", label = "Minimal base mean:", value = 0),
        numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0)
      )
    ),
    mainPanel(
      textOutput("dataAvailabilityMessage"),
      tabsetPanel(
        tabPanel("Summary Table", DTOutput("summaryTable")),
        tabPanel("Manhattan Plot", plotOutput("manhattanPlot")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Dynamically render the ethnicity selection based on analysis type
  output$ethnicity_ui <- renderUI({
    if (input$analysisType == "eQTL") {
      selectInput("ethnicity", "Select Ethnicity:", choices = c("Hispanics", "Non-Hispanic Whites (NHW)"))
    } else if (input$analysisType == "Differential Expression Analysis") {
      selectInput("ethnicity", "Select Ethnicity:", choices = c("MU-BRAIN", "Hispanics", "Non-Hispanic Whites (NHW)"))
    }
  })
  
  # Render UI elements dynamically based on selection
  output$eqtlFilters <- renderUI({
    req(input$ethnicity, input$analysisType)
    
    tagList(

      if (input$eqtlType == "trans") {
        sliderInput("pvalue", "P-Value Threshold:", min = 0, max = 1, value = 0.05, step = 0.01)
      } else if (input$eqtlType == "cis") {
        sliderInput("FDR", "Adjusted P-Value Threshold:", min = 0, max = 1, value = 0.05, step = 0.01)
      },
      actionButton("filterBtn", "Apply Filter")
    )
  })
  
  # Reactive variables for data loading
  eqtl_data <- reactive({
    req(input$analysisType == "eQTL", input$ethnicity, input$eqtlType)
    
    # Debugging messages
    message("eqtl_data() called with:")
    message("Ethnicity: ", input$ethnicity)
    message("eQTL Type: ", input$eqtlType)
    
    if (is.null(input$ethnicity) || is.null(input$eqtlType)) {
      stop("Error: Ethnicity or eQTL Type is NULL in eqtl_data() reactive function")
    }
    
    if (input$eqtlType == "cis"){
      
      withProgress(message = "Loading cis-eQTL data...", value = 0, {
        data <- load_eqtl(package_extdata_path, input$ethnicity, input$eqtlType)
      })
      
    } else {
      
      withProgress(message = "Loading trans-eQTL data...", value = 0, {
        data <- load_eqtl(package_extdata_path, input$ethnicity, input$eqtlType)
      })
      
    }
    
    
    # data <- load_eqtl(package_extdata_path, input$ethnicity, input$eqtlType)
    
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }

  })
  
  
  deseq2_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
    withProgress(message = "Loading DESeq2 data...", value = 0, {
      load_deseq2(package_extdata_path, input$ethnicity)
    })
  })
  
  # Observing and updating gene dropdown for eQTL analysis
  eqtl_data <- reactive({
    req(input$analysisType == "eQTL", input$ethnicity, input$eqtlType)
    
    message("eqtl_data() called with:")
    message("Ethnicity: ", input$ethnicity)
    message("eQTL Type: ", input$eqtlType)
    
    if (is.null(input$ethnicity) || is.null(input$eqtlType)) {
      stop("Error: Ethnicity or eQTL Type is NULL in eqtl_data() reactive function")
    }
    
    data <- load_eqtl(package_extdata_path, input$ethnicity, input$eqtlType)
    
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }
    
    return(data)
  })
  
  # Observing and updating gene dropdown for eQTL analysis
  observe({
    if (input$analysisType == "eQTL") {
      
      if (input$eqtlType == "cis"){
        req(eqtl_data(), input$eqtlType, input$ethnicity, input$FDR)
      } else {
        
        req(eqtl_data(), input$eqtlType, input$ethnicity, input$pvalue)
      }
      
      
      data <- eqtl_data()
      
      if (is.null(data)) {
        updateSelectizeInput(session, "gene", choices = NULL, server = TRUE)
        return()
      }
      
      if (nrow(data) == 0) {
        message("No genes satisfy the filter criteria.")
        updateSelectizeInput(session, "gene", choices = NULL, server = TRUE)
        return()
      }
      
      filtered_genes <- tryCatch({
        if (input$eqtlType == "trans") {
          data %>% filter(pvalue < input$pvalue) %>% pull(gene) %>% unique()
        } else {
          data %>% filter(FDR < input$FDR) %>% pull(gene) %>% unique()
        }
      }, error = function(e) {
        NULL
      })
      
      message("Total = ", length(filtered_genes), " are added in the drop-down list.")
      
      updateSelectizeInput(session, "gene", choices = filtered_genes, server = TRUE)
    } else {
      updateSelectizeInput(session, "gene", choices = NULL, server = TRUE)
    }
  })
  
  # Filtering data based on analysis type
  filtered_data <- reactive({
    if (input$analysisType == "eQTL") {
      
      if (input$eqtlType == "cis"){
        req(eqtl_data(), input$FDR, input$gene)
      } else {
        
        req(eqtl_data(), input$pvalue, input$gene)
      }
      
      message("Filtering data for gene = ", input$gene)
      if (input$eqtlType == "trans") {
        
        eqtl_data() %>% filter(gene == input$gene & pvalue <= input$pvalue)
        
      } else {
        
        eqtl_data() %>% filter(gene == input$gene & FDR <= input$FDR)
        
      }
      
    } else {
      req(input$p_adjusted, input$base_mean, input$log2fc)
      deseq2_data() %>% filter(padj <= as.numeric(input$p_adjusted) & baseMean >= input$base_mean & abs(log2FoldChange) >= input$log2fc)
    }
  })
  
  # Summary table
  output$summaryTable <- renderDT({
    req(filtered_data())
    message("Showing ", dim(filtered_data())[1], " entries for the filtered data....")
    datatable(filtered_data())
  })
  
  
  # Volcano Plot
  output$volcanoPlot <- renderPlot({
    req(input$analysisType == "Differential Expression Analysis", filtered_data())
    res <- filtered_data()
    res$diffexpressed <- ifelse(res$log2FoldChange > 0, "UP", "DOWN")
    colnames(res)[1] <- "gene_symbol"
    res$delabel <- ifelse(res$gene_symbol %in% head(res[order(res$padj), "gene_symbol"], 30), res$gene_symbol, NA)
    
    if  (input$ethnicity == "MU-BRAIN"){
      
      ggplot(data = res, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
        geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
        geom_point(size = 2) +
        # Updated colors for UP, DOWN, and NO genes
        scale_color_manual(values = c("DOWN" = "#00AFBB", "UP" = "#bb0c00"), 
                           labels = c("Downregulated",  "Upregulated")) + 
        coord_cartesian(ylim = c(0, 25), xlim = c(-2, 2)) + 
        labs(color = 'Expression Status', # Updated legend title
             x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
        scale_x_continuous(breaks = seq(-2, 2, 0.5)) + 
        ggtitle(paste0('Volcano Plot for ', input$ethnicity, ' ethnicity')) +
        geom_text_repel(max.overlaps = 20, size = 5, box.padding = 0.5, max.time = 3) + 
        theme_minimal(base_size = 14)
      
    } else {
      
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
        geom_text_repel(max.overlaps = 20, size = 5, box.padding = 0.5, max.time = 3) + 
        theme_minimal(base_size = 14)
      
    }
    
  }, height = 800, width = 800)
  
  
  # Manhattan Plot
  output$manhattanPlot <- renderPlot({
    req(filtered_data())
    
    if (input$analysisType == "eQTL") {
      
      # Manhattan plot for eQTL
      if (input$eqtlType == "trans"){
      filtered_data2 <- filtered_data() %>% rename(FDR = pvalue)
      } else {
        
        filtered_data2 <- filtered_data()
      }
      
      manhattan_df <- data.frame(
        SNP = filtered_data2$snp,
        Chromosome = filtered_data2$CHR,
        Position = filtered_data2$POS,
        P = filtered_data2$FDR
      )
      
      
      if (nrow(manhattan_df) <= 2) {
        output$manhattanPlotMessage <- renderText("The Manhattan plot can't be generated for less than two SNPs.")
        plot(1, type = "n", xlab = "", ylab = "", main = "The Manhattan plot can't be generated for less than two SNPs.")
        return()
      }
      
      output$manhattanPlotMessage <- renderText("")  # Clear any previous message
      
      manhattan_df <- manhattan_df %>% arrange(P)
      if (dim(manhattan_df)[1] > 5){
        
        snp_annot <- manhattan_df$SNP[1:5]
        
      } else {
        
        snp_annot <- manhattan_df$SNP
        
      }
      
      manhattan_df <- manhattan_df %>% arrange(Chromosome, Position)
      
      unique_chr_list <- unique(manhattan_df$Chromosome)
      # message("Setting ", unique_chr_list, " as x-axis labels....")
      CMplot(manhattan_df,
             type="p",
             plot.type="m",
             col = c("blue", "red"),
             LOG10=TRUE,
             file = NULL,
             cex = 0.6, main = paste0("Manhattan Plot of ", input$eqtlType, "-eQTL Analysis for ", input$gene, " Gene in ", input$ethnicity, " Cohort"),
             file.output=F,verbose=F,
             width=14,height=6,chr.labels.angle=45,
             highlight=snp_annot, 
             highlight.text=snp_annot, 
             highlight.text.cex=1.4,
             highlight.col = "red", 
             chr.labels = unique_chr_list,
             chr.den.col = NULL
      )
    }
  })
  
  
}

shinyApp(ui = ui, server = server)
