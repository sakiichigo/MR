# 加载所需的包
library(TwoSampleMR)
library(stringr)
library(dplyr)
library(parallel)

# 读取 eQTL 数据
eqtl_data <- read.table("E:/eqtl/file/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt", header = TRUE)

# 获取 GWAS 汇总数据目录下的所有 txt 文件
gwas_dir <- "D:/document/bioInfo/heart-sarco-0330/outcome"
gwas_files <- list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE)

# 读取 all_common_genes.txt 文件
common_genes_file <- "E:/anaFiles/smrBrainEye/cataractGlaucoma/all_common_genes.txt"
common_genes <- readLines(common_genes_file)

# 数据预处理：筛选显著的 SNP 作为工具变量
eqtl_significant <- eqtl_data[eqtl_data$Pvalue < 5e-8, ]

# 从 eqtl_significant 的 GeneSymbol 列筛选出与 common_genes 重合的内容
eqtl_significant <- eqtl_significant[eqtl_significant$GeneSymbol %in% common_genes, ]

# 根据基因分类拆分 exposure_dat
exposure_dat_list <- split(eqtl_significant, eqtl_significant$GeneSymbol)

exposure_dir <- "E:/eqtl"
if (!dir.exists(exposure_dir)) {
  dir.create(exposure_dir, recursive = TRUE)
}

# 定义处理单个基因数据的函数
process_gene <- function(i, exposure_dat_list, exposure_dir) {
  gene_name <- names(exposure_dat_list[i])
  exposure_file <- file.path(exposure_dir, paste0("exposure_dat_", gene_name, ".csv"))
  
  # 检查文件是否已存在
  if (file.exists(exposure_file)) {
    message(paste("File", exposure_file, "already exists. Skipping..."))
    return(NULL)
  }
  
  exp_dat <- data.frame(
    effect_allele.exposure = exposure_dat_list[[i]]$AssessedAllele,
    other_allele.exposure = exposure_dat_list[[i]]$OtherAllele,
    beta.exposure = exposure_dat_list[[i]]$Zscore,
    se.exposure = 1,  # 当使用 Z-score 时，标准误设为 1
    pval.exposure = exposure_dat_list[[i]]$Pvalue,
    SNP = exposure_dat_list[[i]]$SNP,
    exposure = exposure_dat_list[[i]]$GeneSymbol,
    id.exposure = exposure_dat_list[[i]]$Gene
  )
  
  # LD 聚类
  exposure_dat_clumped <- tryCatch({
    ld_clump(
      dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure),
      plink_bin = "D:/program/plink/plink.exe",
      bfile = "D:/document/bioInfo/EUR/EUR"
    )
  }, error = function(e) {
    print(e)
    return(NULL)
  })
  
  if (is.null(exposure_dat_clumped)) {
    return(NULL)
  }
  
  exposure_dat <- exp_dat[exp_dat$SNP %in% exposure_dat_clumped$rsid, ]
  
  # 死循环
  while (TRUE) {
    error_occurred <- FALSE
    
    if (is.null(exposure_dat$eaf.exposure) || any(is.na(exposure_dat$eaf.exposure))) {
      exposure_dat_add <- tryCatch({
        snp_add_eaf(exposure_dat)
      }, error = function(e) {
        # 检查是否为特定错误信息
        error_messages <- c(
          "seaching for No.\\d+ SNP",
          "Bad Gateway \\(HTTP 502\\)",
          "列表参数有错误：所有变量的长度都应该是一样的"
        )
        
        if (any(sapply(error_messages, grepl, toString(e)))) {
          cat("捕获到指定错误，继续执行...\n")
        } else {
          print(e)
        }
        
        error_occurred <<- TRUE
        return(NULL)
      })
    }
    
    # 如果没有错误发生，跳出循环
    if (!error_occurred) {
      break
    }
  }
  
  # 写入文件
  if (!is.null(exposure_dat_add)) {
    write.csv(exposure_dat_add, exposure_file, row.names = FALSE)
  } else {
    write.csv(exposure_dat, exposure_file, row.names = FALSE)
  }
  
  return(NULL)
}

# 获取可用的核心数
num_cores <- detectCores()

# 创建并行集群
cl <- makeCluster(num_cores)

# 向集群节点导出必要的变量和函数
clusterExport(cl, c("exposure_dat_list", "exposure_dir", "ld_clump", "snp_add_eaf", "tibble"))

# 并行处理每个基因数据
parLapply(cl, seq_along(exposure_dat_list), process_gene, exposure_dat_list, exposure_dir)

# 停止集群
stopCluster(cl)

# 后续代码保持不变
# 初始化错误数据向量
errorData_out <- c()
errorData_exp <- c()

# 结果保存目录
result_dir <- "D:/document/bioInfo/heart-sarco-0330/eqtl/result"
# 创建结果目录（如果不存在）
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

# 获取 E:/eqtl 目录下的所有 csv 文件
exposure_files <- list.files(exposure_dir, pattern = "\\.csv$", full.names = TRUE)
FLAG=TRUE
#并行
if(FALSE){
  
  # 提前读取所有的 GWAS 文件
  gwas_data_list <- lapply(gwas_files, function(gwas_file) {
    read.table(gwas_file, header = FALSE, sep = "\t", fill = TRUE)
  })
  
  # 并行处理函数
  process_exposure_file <- function(exposure_file, exposure_index, total_exposure_files) {
    cat(paste0("开始处理第 ", exposure_index, " 个暴露文件（共 ", total_exposure_files, " 个）：", basename(exposure_file), "\n"))
    exp_dat <- read.csv(exposure_file)
    if (length(exp_dat) == 0) {
      cat(paste0("第 ", exposure_index, " 个暴露文件为空，跳过。\n"))
      return(NULL)
    }
    gene_name <- sub("^.*_", "", sub("\\.csv$", "", basename(exposure_file)))
    
    for (i in seq_along(gwas_files)) {
      gwas_data <- gwas_data_list[[i]]
      unique_ids <- unique(gwas_data$V1)
      total_ids <- length(unique_ids)
      for (id_index in seq_along(unique_ids)) {
        id <- unique_ids[id_index]
        cat(paste0("  处理第 ", i, " 个 GWAS 文件中的第 ", id_index, " 个 ID（共 ", total_ids, " 个）：", id, "\n"))
        # 构造结果文件路径
        dat_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_dat.csv"))
        res_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_res.csv"))
        or_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_or.csv"))
        
        # 检查文件是否存在，如果存在则跳过本次循环
        if (file.exists(dat_file) && file.exists(res_file) && file.exists(or_file)) {
          cat(paste0("    结果文件已存在，跳过。\n"))
          next
        }
        
        # 初始化标志变量
        message_to_next <- TRUE
        error_to_next <- FALSE
        
        while (TRUE) {
          result <- tryCatch({
            withCallingHandlers(
              outcome_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = id),
              message = function(c) if (stringr::str_detect(as.character(c), "Failed to"))
                message_to_next <<- FALSE
            )
            error_to_next <- TRUE
          },
          error = function(e) {
            print(e)
            if (grepl("used up your OpenGWAS allowance", as.character(e), fixed = TRUE)) {
              message("检测到达到 OpenGWAS 配额限制，等待 5 分钟后继续...")
              Sys.sleep(300) # 等待 5 分钟（300 秒）
            }
            return(NULL)
          })
          
          if (message_to_next == TRUE & error_to_next == TRUE) {
            break
          }
          if (is.null(result)) {
            outcome_dat <- NULL
            break
          }
        }
        
        # 数据匹配
        if (!is.null(outcome_dat)) {
          dat <- harmonise_data(exp_dat, outcome_dat)
          
          if (!length(dat) || dim(dat)[1] == 0) {
            errorData_out <<- c(errorData_out, id)
            errorData_exp <<- c(errorData_exp, exp_dat$SNP)
            cat(paste0("    数据匹配后无有效数据，跳过。\n"))
            next
          }
          
          # 异质性多效性
          ple <- mr_pleiotropy_test(dat)
          he <- mr_heterogeneity(dat)
          
          if (length(ple$pval) == 0 || length(he$Q_pval) == 0) {
            # 可能出现没有足够 SNP
            res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
            ivw <- "Inverse variance weighted"
          } else {
            if (he[he$method == "Inverse variance weighted", ]$Q_pval < 0.05) {
              he <- mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre")) # he 更改方法
              # 存在异质性时随机效应 ivw_mre，同时 raps 替换 simple mode
              res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_raps_modified", "mr_weighted_mode"))
              ivw <- "Inverse variance weighted (multiplicative random effects)"
            } else {
              res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
              ivw <- "Inverse variance weighted"
            }
          }
          
          # 空数据
          if (nrow(res) == 0) {
            cat(paste0("    分析结果为空，跳过。\n"))
            next
          }
          
          res[is.na(res$method), ]$method <- "Robust Adjusted Profile Score"
          or <- generate_odds_ratios(res)
          
          # or 写入 HE，只保留 IVW
          if (length(he) == 0) {
            or[, 15] <- NA
          } else {
            if (dim(he)[1] > 1 & he[2, ]$method == ivw) {
              or[, 15] <- he[he$method == ivw, ]$Q_pval
            } else {
              or[, 15] <- he[1, ]$Q_pval
            }
          }
          colnames(or)[15] <- "he"
          
          # or 写入 PLE
          if (length(ple) == 0) {
            or[, 16] <- NA
          } else {
            or[, 16] <- ple[1, 7]
          }
          colnames(or)[16] <- "ple"
          
          or[, 17] <- Ff(exp_dat)$fm
          colnames(or)[17] <- "F"
          
          write.csv(dat, dat_file, row.names = FALSE)
          write.csv(res, res_file, row.names = FALSE)
          write.csv(or, or_file, row.names = FALSE)
          cat(paste0("    结果文件已保存。\n"))
        }
      }
    }
    cat(paste0("第 ", exposure_index, " 个暴露文件处理完成。\n"))
  }
  
  # 并行处理所有的暴露文件
  num_cores <- detectCores() - 1  # 使用除一个核心外的所有核心
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, {
    library(TwoSampleMR)
    library(stringr)
    sink(stdout())  # 将输出重定向到主会话
  })
  clusterExport(cl, c("exposure_files", "gwas_files", "gwas_data_list", "result_dir", "errorData_out", "errorData_exp", "Ff", "extract_outcome_data", "harmonise_data", "mr_pleiotropy_test", "mr_heterogeneity", "mr", "generate_odds_ratios"))
  
  total_exposure_files <- length(exposure_files)
  for (i in seq_along(exposure_files)) {
    exposure_file <- exposure_files[i]
    parLapply(cl, list(exposure_file), process_exposure_file, exposure_index = i, total_exposure_files = total_exposure_files)
  }
  
  stopCluster(cl)
}else{
  
  # 对 E:/eqtl 目录下的 csv 文件进行循环
  for (exposure_file in exposure_files) {
    
    exp_dat <- read.csv(exposure_file)
    if (length(exp_dat) == 0) {
      next  # 跳过本次循环，进入下一次循环
    }
    gene_name <- sub("^.*_", "", sub("\\.csv$", "", basename(exposure_file)))
    
    # 读取当前的 GWAS 汇总数据
    for (gwas_file in gwas_files) {
      gwas_data <- read.table(gwas_file, header = FALSE, sep = "\t")
      
      # 假设 gwas_data 中有 id 列，对每个 id 进行循环
      for (id in unique(gwas_data$V1)) {
        # 构造结果文件路径
        dat_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_dat.csv"))
        res_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_res.csv"))
        or_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_or.csv"))
        
        # 检查文件是否存在，如果存在则跳过本次循环
        if (file.exists(dat_file) && file.exists(res_file) && file.exists(or_file)) {
          next
        }
        
        # 初始化标志变量
        message_to_next <- TRUE
        error_to_next <- FALSE
        
        while (TRUE) {
          result <- tryCatch({
            withCallingHandlers(
              outcome_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = id),
              message = function(c) if (stringr::str_detect(as.character(c), "Failed to"))
                message_to_next <<- FALSE
            )
            error_to_next <- TRUE
          },
          error = function(e) {
            print(e)
            if (grepl("used up your OpenGWAS allowance", as.character(e), fixed = TRUE)) {
              message("检测到达到 OpenGWAS 配额限制，等待 5 分钟后继续...")
              Sys.sleep(300) # 等待 5 分钟（300 秒）
            }
            return(NULL)
          })
          
          if (message_to_next == TRUE & error_to_next == TRUE) {
            break
          }
          if (is.null(result)) {
            outcome_dat <- NULL
            break
          }
        }
        
        # 数据匹配
        if (!is.null(outcome_dat)) {
          dat <- harmonise_data(exp_dat, outcome_dat)
          
          if (!length(dat) || dim(dat)[1] == 0) {
            errorData_out <- c(errorData_out, id)
            errorData_exp <- c(errorData_exp, exp_dat$SNP)
            next
          }
          
          # 异质性多效性
          ple <- mr_pleiotropy_test(dat)
          he <- mr_heterogeneity(dat)
          
          if (length(ple$pval) == 0 || length(he$Q_pval) == 0) {
            # 可能出现没有足够 SNP
            res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
            ivw <- "Inverse variance weighted"
          } else {
            if (he[he$method == "Inverse variance weighted", ]$Q_pval < 0.05) {
              he <- mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre")) # he 更改方法
              # 存在异质性时随机效应 ivw_mre，同时 raps 替换 simple mode
              res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_raps_modified", "mr_weighted_mode"))
              ivw <- "Inverse variance weighted (multiplicative random effects)"
            } else {
              res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
              ivw <- "Inverse variance weighted"
            }
          }
          
          # 空数据
          if (nrow(res) == 0) {
            next
          }
          
          res[is.na(res$method), ]$method <- "Robust Adjusted Profile Score"
          or <- generate_odds_ratios(res)
          
          # or 写入 HE，只保留 IVW
          if (length(he) == 0) {
            or[, 15] <- NA
          } else {
            if (dim(he)[1] > 1 & he[2, ]$method == ivw) {
              or[, 15] <- he[he$method == ivw, ]$Q_pval
            } else {
              or[, 15] <- he[1, ]$Q_pval
            }
          }
          colnames(or)[15] <- "he"
          
          # or 写入 PLE
          if (length(ple) == 0) {
            or[, 16] <- NA
          } else {
            or[, 16] <- ple[1, 7]
          }
          colnames(or)[16] <- "ple"
          
          or[, 17] <- Ff(exp_dat)$fm
          colnames(or)[17] <- "F"
          
          write.csv(dat, dat_file, row.names = FALSE)
          write.csv(res, res_file, row.names = FALSE)
          write.csv(or, or_file, row.names = FALSE)
        }
      }
    }
  }
}
