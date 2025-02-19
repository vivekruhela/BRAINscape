#' Load eQTL data
#' @param path The base path to data files
#' @param ethnicity "Hispanics" or "NonHispanics"
#' @param type "cis" or "trans"
#' @return The loaded eQTL data
#' @export
load_eqtl <- function(path, ethnicity, type) {
  
  # Check if ethnicity or type is NULL or missing
  if (missing(ethnicity) || is.null(ethnicity) || ethnicity == "") {
    stop("Error: 'ethnicity' is missing or NULL in load_eqtl()")
  }
  
  if (missing(type) || is.null(type) || type == "") {
    stop("Error: 'type' is missing or NULL in load_eqtl()")
  }
  
  message("Loading eQTL data for ethnicity: ", ethnicity, " and type: ", type)
  
  if (type == "cis") {
    if (ethnicity == "Hispanics") {
      tryCatch({
        data <- readRDS(file.path(path, "eqtl/Hispanics/BRAIN_Hispanics_tensoreqtl_cis_eqtl_allchr.rds"))
        message("Data loaded successfully for Hispanics cis-eQTL with ", dim(data)[1], " rows and ", dim(data)[2], " columns.")
      }, error = function(e) { stop("Error loading data: ", e$message) })
      
    } else {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_cis_eqtl_allchr.rds"))
      message("Data loaded successfully for Hispanics trans-eQTL with ", dim(data)[1], " rows and ", dim(data)[2], " columns.")
    }
  } else {  # type == "trans"
    if (ethnicity == "Non-Hispanic Whites (NHW)") {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_trans_eqtl_allchr.rds"))
      message("Data loaded successfully for NHW cis-eQTL with ", dim(data)[1], " rows and ", dim(data)[2], " columns.")
    } else {
      data <- readRDS(file.path(path, "eqtl/Hispanics/BRAIN_Hispanics_tensoreqtl_trans_eqtl_allchr.rds"))
      message("Data loaded successfully for NHW trans-eQTL with ", dim(data)[1], " rows and ", dim(data)[2], " columns.")
    }
  }
  
  # Ensure data is not NULL
  if (is.null(data)) {
    stop("Error: eQTL data could not be loaded. Check file paths.")
  }
  
  # Filter out unwanted genes
  data <- data %>% filter(!grepl("^MIR", gene)) %>% filter(!grepl("^LINC", gene))
  
  data
}

#' Load DESeq2 data
#' @export
load_deseq2 <- function(path, ethnicity) {
  if (ethnicity == "Hispanics") {
    data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_hispanic.rds"))
  } else {
    if (ethnicity == "Non-Hispanic Whites (NHW)") {
      data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_white.rds"))
    } else {
      data <- readRDS(file.path(path, "deseq2/DESeq2_model1_results.rds"))
    }
    
  }
    data
}
