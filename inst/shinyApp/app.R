# Load required libraries
suppressPackageStartupMessages({
  
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(CMplot)
  library(DT)
  library(ggrepel)
  library(plotly)
  
})

rm(list = ls())
gc()

# Load functions
package_extdata_path <- paste(system.file(package = "BRAINscape"), "extdata", sep = "/")

ui <- fluidPage(
  titlePanel("Analysis Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("analysisType", "Select Analysis Type:", 
                  choices = c("eQTL", "Differential Expression Analysis", "GWAS")),
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
        radioButtons("BRAAK", "Select BRAAK Type:", choices = c("Binary", "Continuous")),
        radioButtons("Model", "Select Model:", choices = c("Model1", "Model2", "Model3", "Model4")),
        selectInput("p_adjusted", label = "Adjusted P-Value Threshold:", 
                    choices = c("0.00001" = 0.00001,"0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0)),
        numericInput("base_mean", label = "Minimal base mean:", value = 0),
        numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0)
      ),
      
      # Conditional UI for GWAS analysis
      conditionalPanel(
        condition = "input.analysisType == 'GWAS'",
        selectInput("P", label = "P-Value Threshold:", 
                    choices = c("0.00001" = 0.00001,"0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0))
      )
    ),
    mainPanel(
      textOutput("dataAvailabilityMessage"),
      tabsetPanel(
        tabPanel("Summary Table", DTOutput("summaryTable")),
        tabPanel("Manhattan Plot", plotOutput("manhattanPlot")),
        tabPanel("Volcano Plot", plotlyOutput("volcanoPlot"))
        
      )
    )
  )
)

server <- function(input, output, session) {
  # Dynamically render the ethnicity selection based on analysis type
  output$ethnicity_ui <- renderUI({
    if (input$analysisType == "eQTL" || input$analysisType == "GWAS") {
      selectInput("ethnicity", "Select Ethnicity:", choices = c("Hispanics", "Non-Hispanic Whites (NHW)"))
    } else if (input$analysisType == "Differential Expression Analysis") {
      selectInput("ethnicity", "Select Ethnicity:", choices = c("Hispanics", "Non-Hispanic Whites (NHW)", "MUBRAIN"))
    }
  })
  
  # Render UI elements dynamically based on selection
  output$eqtlFilters <- renderUI({
    req(input$ethnicity, input$analysisType)
    tagList(
      if (input$eqtlType == "trans") {
        selectInput("pval", label = "P-Value Threshold:", choices = c("0.00001" = 0.00001,"0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0))
      } else if (input$eqtlType == "cis") {
        selectInput("qval", label = "Adjusted P-Value Threshold:", choices = c("0.0001" = 0.0001,"0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "1.0" = 1.0))
      }
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
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }
  })
  
  deseq2_data <- reactive({
    req(input$analysisType == "Differential Expression Analysis")
    withProgress(message = "Loading DESeq2 data...", value = 0, {
      load_deseq2(
        path = package_extdata_path,
        ethnicity = input$ethnicity,
        model = input$Model,
        braak = input$BRAAK
      )
    })
  })
  
  gwas_data <- reactive({
    req(input$analysisType == "GWAS")
    withProgress(message = "Loading GWAS data...", value = 0, {
      load_gwas(
        path = package_extdata_path,
        ethnicity = input$ethnicity
      )
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
        req(eqtl_data(), input$eqtlType, input$ethnicity, input$qval)
      } else {
        
        req(eqtl_data(), input$eqtlType, input$ethnicity, input$pval)
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
          data %>% filter(pval < as.numeric(input$pval)) %>% pull(phenotype_id) %>% unique()
        } else {
          data %>% filter(qval < as.numeric(input$qval)) %>% pull(phenotype_id) %>% unique()
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
      if (input$eqtlType == "cis") {
        req(eqtl_data(), input$qval, input$gene)
        eqtl_data() %>% 
          filter(phenotype_id == input$gene & qval <= as.numeric(input$qval))
      } else {
        req(eqtl_data(), input$pval, input$gene)
        eqtl_data() %>% 
          filter(phenotype_id == input$gene & pval <= as.numeric(input$pval))
      }
    } else {
      if (input$analysisType == "Differential Expression Analysis"){
        req(input$p_adjusted, input$base_mean, input$log2fc)
        deseq2_data() %>% 
          filter(padj <= as.numeric(input$p_adjusted) & 
                   baseMean >= input$base_mean & 
                   abs(log2FoldChange) >= input$log2fc)
      } else {
        req(input$P)
        gwas_data() %>% filter(P <= as.numeric(input$P)) %>% arrange(P)
      }
    }
  })
  
  
  # Summary table
  output$summaryTable <- renderDT({
    req(filtered_data())
    message("Showing ", dim(filtered_data())[1], " entries for the filtered data....")
    datatable(filtered_data())
  })
  
  
  # Volcano Plot
  
  output$volcanoPlot <- renderPlotly({
    req(input$analysisType == "Differential Expression Analysis", filtered_data())
    res <- filtered_data()
    
    # Thresholds
    lfc_thresh <- 0.15
    pval_thresh <- 0.05
    
    # Categorize points
    res$category <- "NS"
    res$category[res$padj < pval_thresh & abs(res$log2FoldChange) >= lfc_thresh & res$log2FoldChange > 0] <- "Up"
    res$category[res$padj < pval_thresh & abs(res$log2FoldChange) >= lfc_thresh & res$log2FoldChange < 0] <- "Down"
    res$category[res$padj < pval_thresh & abs(res$log2FoldChange) < lfc_thresh] <- "Significant"
    
    # Dynamic Y-axis based on data
    maxY <- -log10(min(res$padj, na.rm = TRUE)) + 1
    
    # Build ggplot
    p <- ggplot(data = res, aes(
      x = log2FoldChange,
      y = -log10(padj),
      col = category,
      text = paste0(
        "Gene: ", symbol, "<br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "-log10(padj): ", round(-log10(padj), 3), "<br>",
        "padj: ", signif(padj, 3)
      )
    )) +
      geom_point(size = 2) +
      geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), col = "black", linetype = 'dashed') +
      geom_hline(yintercept = -log10(pval_thresh), col = "black", linetype = 'dashed') +
      scale_color_manual(
        values = c(
          "Up" = "#0072B2",
          "Down" = "#D55E00",
          "Significant" = "#009E73",
          "NS" = "#BEBEBE"
        )
      ) +
      coord_cartesian(ylim = c(0, maxY), xlim = c(-2, 2)) +
      labs(
        x = "log2 Fold Change",
        y = "-log10(adjusted p-value)",
        color = "Category"
      ) +
      scale_x_continuous(breaks = seq(-2, 2, 0.5)) +
      ggtitle(paste0('Volcano Plot for ', input$ethnicity, ' ethnicity')) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        text = element_text(size = 14)
      )
    
    # Convert to interactive plot with fixed size
    ggplotly(p, tooltip = "text") %>%
      layout(
        width = 900,
        height = 900,
        title = list(text = paste0('Differential Expression – ', input$ethnicity, ' | ', input$Model, ' | ', input$BRAAK, ' BRAAK')),
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(p-adjusted value)")
      )
  })
  
  
  # Manhattan Plot
  output$manhattanPlot <- renderPlot({
    req(filtered_data())
    if (input$analysisType == "eQTL") {
      filtered_data2 <- if (input$eqtlType == "trans") {
        filtered_data() %>% rename(FDR = pval)
      } else {
        filtered_data2 <- filtered_data() %>% rename(FDR = qval)
      }
      
      # Build dataframe
      manhattan_df <- data.frame(
        SNP = filtered_data2$variant_id,
        Chromosome = as.factor(filtered_data2$chr),
        Position = filtered_data2$pos,
        P = filtered_data2$FDR
      )
      
      if (nrow(manhattan_df) <= 2) {
        plot.new()
        title("Not enough SNPs to plot.")
        return()
      }
      
      # Compute -log10(p)
      manhattan_df$logP <- -log10(as.numeric(manhattan_df$P))
      manhattan_df <- manhattan_df %>% arrange(as.numeric(Chromosome), Position)
      chrom_sizes <- manhattan_df %>% 
        group_by(Chromosome) %>%
        summarise(chr_len = max(Position), .groups = "drop")
      
      chrom_sizes$cum_len <- cumsum(chrom_sizes$chr_len) - chrom_sizes$chr_len
      manhattan_df <- manhattan_df %>%
        left_join(chrom_sizes %>% select(Chromosome, cum_len), by = "Chromosome") %>%
        mutate(Pos_cum = Position + cum_len)
      
      # Axis ticks (midpoints per chromosome)
      axis_df <- manhattan_df %>%
        group_by(Chromosome) %>%
        summarise(center = mean(range(Pos_cum)), .groups = "drop")
      
      # Select top SNPs for annotation
      top_snps <- manhattan_df %>% arrange(P) %>% slice_head(n = 5)
      
      ggplot(manhattan_df, aes(x = Pos_cum, y = logP, color = as.factor(as.numeric(Chromosome) %% 2))) +
        geom_point(alpha = 0.8, size = 1.5) +
        scale_color_manual(values = c("#2111a1", "#a11111"), guide = "none") +
        geom_point(data = top_snps, aes(x = Pos_cum, y = logP), color = "red", size = 3) +
        ggrepel::geom_text_repel(
          data = top_snps,
          aes(x = Pos_cum, y = logP, label = SNP),
          size = 4,
          nudge_y = 0.5,
          segment.color = "black"
        ) +
        scale_x_continuous(
          label = axis_df$Chromosome,
          breaks = axis_df$center
        ) +
        labs(
          title = paste0("Manhattan Plot: ", input$phenotype_id, " in ", input$ethnicity),
          x = "Chromosome",
          y = "-log10(p)"
        ) +
        theme_minimal(base_size = 16) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
    } else {
      message("Generating GWAS Manhattan Plot for ", dim(filtered_data())[1], " SNPs.")
      
      pval_thresh <- 5E-08
      
      # Build dataframe
      manhattan_df <- data.frame(
        SNP = filtered_data()$ID,
        Chromosome = filtered_data()$CHROM,
        Position = filtered_data()$POS,
        P = filtered_data()$P
      )
      
      # Compute -log10(p)
      manhattan_df$P[manhattan_df$P == 0] <- 1e-300
      manhattan_df$logP <- -log10(as.numeric(manhattan_df$P))
      
      # Select top SNPs for annotation
      top_snps <- manhattan_df %>% arrange(P) %>% slice_head(n = 5)
      
      manhattan_df <- manhattan_df %>% arrange(as.numeric(Chromosome), Position)
      chrom_sizes <- manhattan_df %>% 
        group_by(Chromosome) %>%
        summarise(chr_len = max(Position), .groups = "drop")
      
      # chrom_sizes$cum_len <- cumsum(chrom_sizes$chr_len) - chrom_sizes$chr_len
      chrom_sizes$cum_len <- cumsum(as.numeric(chrom_sizes$chr_len)) - as.numeric(chrom_sizes$chr_len)
      
      manhattan_df <- manhattan_df %>%
        left_join(chrom_sizes %>% select(Chromosome, cum_len), by = "Chromosome") %>%
        mutate(Pos_cum = Position + cum_len)
      
      # Axis ticks (midpoints per chromosome)
      axis_df <- manhattan_df %>%
        group_by(Chromosome) %>%
        summarise(center = mean(range(Pos_cum)), .groups = "drop")
      
      # Select top SNPs for annotation
      top_snps <- manhattan_df %>% arrange(P) %>% slice_head(n = 5)
      
      ggplot(manhattan_df, aes(x = Pos_cum, y = logP, color = as.factor(as.numeric(Chromosome) %% 2))) +
        geom_point(alpha = 0.8, size = 1.5) +
        scale_color_manual(values = c("#2111a1", "#a11111"), guide = "none") +
        geom_point(data = top_snps, aes(x = Pos_cum, y = logP), color = "red", size = 3) +
        ggrepel::geom_text_repel(
          data = top_snps,
          aes(x = Pos_cum, y = logP, label = SNP),
          size = 4,
          nudge_y = 0.5,
          segment.color = "black"
        ) +
        geom_hline(yintercept = -log10(pval_thresh), col = "#0c7326", linetype = 'dashed') +
        scale_x_continuous(
          label = axis_df$Chromosome,
          breaks = axis_df$center
        ) +
        labs(
          title = paste0("Manhattan Plot: ", input$ethnicity, " (p-value < ", input$P, ")"),
          x = "Chromosome",
          y = "-log10(p)"
        ) +
        theme_minimal(base_size = 16) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }
  })
  
  
}

shinyApp(ui = ui, server = server)
