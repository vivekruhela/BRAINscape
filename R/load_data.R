#' Load eQTL data
#' @param path The base path to data files
#' @param ethnicity "Hispanics" or "NonHispanics" or "MUBRAIN"
#' @param type "cis" or "trans" or "Model"
#' @return The loaded data
#' @return A data frame containing GWAS results.
#' @export
load_eqtl <- function(path, ethnicity, type) {
  
  if (ethnicity == "Hispanics") {
    if (type == "cis"){
      data <- readRDS(file.path(path, "eqtl", "Hispanics/BRAIN.cis_qtl_results.merged.rds"))
    } else {
      data <- readRDS(file.path(path, "eqtl","Hispanics/BRAIN.trans_qtl_results.merged_top1M.rds"))
    }
  }
  
  if (ethnicity == "Non-Hispanic Whites (NHW)") {
    if (type == "cis"){
      data <- readRDS(file.path(path, "eqtl", "NonHispanics/BRAIN.cis_qtl_results.500kb_0.05.merged2.rds"))
    } else {
      data <- readRDS(file.path(path, "eqtl", "NonHispanics/BRAIN.trans_qtl_results.500kb_0.05.merged_top1M.rds"))
    }
  }
  
  message(type, "-eQTL analysis data for ", ethnicity, "is  loaded successfully with columns = ", paste(colnames(data), " "), " and rows = ", nrow(data), ".")
  data
}


#' Load DESeq2 data
#' @export
load_deseq2 <- function(path, ethnicity, model, braak) {
  
  if (ethnicity == "Hispanics") {
    if (model == 'Model1') {
      if (braak == 'Binary') {
        data <- readRDS(file.path(path, "deseq2/Model1/model1_N113_Hispanics_only_binary_braak.rds"))
      } else {
        data <- readRDS(file.path(path, "deseq2/Model1/model1_N113_Hispanics_only_conti_braak.rds"))
      }
    }
    
    if (model == 'Model2') {
      if (braak == 'Binary') {
        data <- readRDS(file.path(path, "deseq2/Model2/model2_N69_Hispanics_only_binary_braak.rds"))
      } else {
        data <- readRDS(file.path(path, "deseq2/Model2/model2_N69_Hispanics_only_conti_braak.rds"))
      }
    }
    
    if (model == 'Model3') {
      if (braak == 'Binary') {
        data <- readRDS(file.path(path, "deseq2/Model3/model3_N136_Hispanics_only_binary_braak.rds"))
      } else {
        data <- readRDS(file.path(path, "deseq2/Model3/model3_N136_Hispanics_only_conti_braak.rds"))
      }
    }
    
    if (model == 'Model4') {
      if (braak == 'Binary') {
        data <- readRDS(file.path(path, "deseq2/Model4/model4_N89_Hispanics_only_binary_braak.rds"))
      } else {
        data <- readRDS(file.path(path, "deseq2/Model4/model4_N89_Hispanics_only_conti_braak.rds"))
      }
    }
    
  } else {
    if (ethnicity == "Non-Hispanic Whites (NHW)"){
      
      if (model == 'Model1') {
        if (braak == 'Binary') {
          data <- readRDS(file.path(path, "deseq2/Model1/model1_N399_NHW_only_binary_braak.rds"))
        } else {
          data <- readRDS(file.path(path, "deseq2/Model1/model1_N399_NHW_only_conti_braak.rds"))
        }
      }
      
      if (model == 'Model2') {
        if (braak == 'Binary') {
          data <- readRDS(file.path(path, "deseq2/Model2/model2_N203_NHW_only_binary_braak.rds"))
        } else {
          data <- readRDS(file.path(path, "deseq2/Model2/model2_N203_NHW_only_conti_braak.rds"))
        }
      }
      
      if (model == 'Model3') {
        if (braak == 'Binary') {
          data <- readRDS(file.path(path, "deseq2/Model3/model3_N512_NHW_only_binary_braak.rds"))
        } else {
          data <- readRDS(file.path(path, "deseq2/Model3/model3_N512_NHW_only_conti_braak.rds"))
        }
      }
      
      if (model == 'Model4') {
        if (braak == 'Binary') {
          data <- readRDS(file.path(path, "deseq2/Model4/model4_N309_NHW_only_binary_braak.rds"))
        } else {
          data <- readRDS(file.path(path, "deseq2/Model4/model4_N309_NHW_only_conti_braak.rds"))
        }
      }
    } else {
      if (ethnicity == "MUBRAIN"){
        
        if (model == 'Model1') {
          if (braak == 'Binary') {
            data <- readRDS(file.path(path, "deseq2/Model1/model1_N565_binary_braak.rds"))
          } else {
            data <- readRDS(file.path(path, "deseq2/Model1/model1_N565_conti_braak.rds"))
          }
        }
        
        if (model == 'Model2') {
          if (braak == 'Binary') {
            data <- readRDS(file.path(path, "deseq2/Model2/model2_N272_binary_braak.rds"))
          } else {
            data <- readRDS(file.path(path, "deseq2/Model2/model2_N272_conti_braak.rds"))
          }
        }
        
        if (model == 'Model3') {
          if (braak == 'Binary') {
            data <- readRDS(file.path(path, "deseq2/Model3/model3_N720_binary_braak.rds"))
          } else {
            data <- readRDS(file.path(path, "deseq2/Model3/model3_N720_conti_braak.rds"))
          }
        }
        
        if (model == 'Model4') {
          if (braak == 'Binary') {
            data <- readRDS(file.path(path, "deseq2/Model4/model4_N398_binary_braak.rds"))
          } else {
            data <- readRDS(file.path(path, "deseq2/Model4/model4_N398_conti_braak.rds"))
          }
        }
      }
    }
  }
  data <- data %>% filter(!is.na(padj) & padj > 0) %>% arrange(padj)
  message("DESeq2 for ", ethnicity, " for ", braak, " BRAAK and model ", model, " results are loaded data successfully with columns = ", paste(colnames(data), " "), " and rows = ", nrow(data), ".")
  data
}

#' Load GWAS data
#' @export
load_gwas <- function(path, ethnicity) {
  
  if (ethnicity == "Hispanics") {
    data <- readRDS(file.path(path, "gwas", "Hispanics/hispancis.BRAAK_presidual.glm.linear.ADD_summary_stats.rds"))
  }
  
  if (ethnicity == "Non-Hispanic Whites (NHW)") {
    data <- readRDS(file.path(path, "gwas", "NonHispanics/nhw.BRAAK_presidual.glm.linear.ADD_summary_stats.rds"))
  }
  
  message("GWAS analysis data for ", ethnicity, "is  loaded successfully with columns = ", paste(colnames(data), " "), " and rows = ", nrow(data), ".")
  data
}

