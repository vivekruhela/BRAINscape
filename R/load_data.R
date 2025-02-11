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
        message("Data loaded successfully")
      }, error = function(e) { message("Error loading data:", e$message) }
      )
      
    } else {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_cis_eqtl_allchr.rds"))
    }
  } else {
    if (ethnicity == "NonHispanics") {
      data <- readRDS(file.path(path, "eqtl/NonHispanics/BRAIN_NonHispanics_tensoreqtl_trans_eqtl_allchr.rds"))
    } else {
      data <- readRDS(file.path(path, "eqtl/Hispanics/BRAIN_Hispanics_tensoreqtl_trans_eqtl_allchr.rds"))
    }
  }
  data
}

#' Load GWAS data
#' @export
load_gwas <- function(path, ethnicity) {
  if (ethnicity == "Hispanics") {
    data <- readRDS(file.path(path, "gwas/Hispanics_allchr_logistic.BRAAK_presidual.glm.linear4.rds"))
  } else {
    data <- readRDS(file.path(path, "gwas/NonHispanics_allchr_logistic.BRAAK_presidual.glm.linear4.rds"))
  }
  data
}
  
#' Load DESeq2 data
#' @export
load_deseq2 <- function(path, ethnicity) {
  if (ethnicity == "Hispanics") {
    data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_hispanic.rds"))
  } else {
    data <- readRDS(file.path(path, "deseq2/DESeq2_Model1_AD_Conti_white.rds"))
  }
    data
}
