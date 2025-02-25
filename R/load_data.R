#' Load eQTL data
#' @param path The base path to data files
#' @param ethnicity "Hispanics" or "NonHispanics"
#' @param type "cis" or "trans"
#' @return The loaded eQTL data
#' @export
load_eqtl <- function(path, ethnicity, type) {

  if (type == "cis") {
    if (ethnicity == "Hispanics") {
      tryCatch( {

        data <- readRDS(file.path(path, "eqtl/Hispanics/BRAIN_Hispanics_tensoreqtl_cis_eqtl_allchr.rds"))
        message("Hispanics cis-eqtl Data loaded successfully with columns = ", colnames(data))
      }, error = function(e) { message("Error loading data:", e$message) }
      )
      
    } else {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_cis_eqtl_allchr.rds"))
      message("NHW cis-eqtl Data loaded successfully with columns = ", colnames(data))
    }
  } else {
    if (ethnicity == "Non-Hispanic Whites (NHW)") {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_trans_eqtl_allchr.rds"))
      message("NHW tras-eqtl Data loaded successfully with columns = ", colnames(data))
    } else {
      data <- readRDS(file.path(path, "eqtl/Hispanics/BRAIN_Hispanics_tensoreqtl_trans_eqtl_allchr.rds"))
      message("Hispanics cis-eqtl Data loaded successfully with columns = ", colnames(data))
    }
  }
  data
}
  
#' Load DESeq2 data
#' @export
load_deseq2 <- function(path, ethnicity) {
  if (ethnicity == "Hispanics") {
    data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_hispanic.rds"))
    message("DESeq2 for Hispanics Brains results are loaded data successfully.")
  } else {
    if (ethnicity == "Non-Hispanic Whites (NHW)"){

      data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_white.rds"))
      message("DESeq2 for NHW Brains results are loaded data successfully.")
    } else {
      data <- readRDS(file.path(path, "deseq2/DESeq2_model1_results.rds"))
      message("DESeq2 results for MU-Brains are loaded data successfully.")
    }
    
  }
    data
}