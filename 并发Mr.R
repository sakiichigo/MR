# 加载 parallel 包
library(parallel)

# 检测操作系统类型
is_windows <- .Platform$OS.type == "windows"

# 对 E:/eqtl 目录下的 csv 文件进行循环
for (i in seq_along(exposure_files)) {
  exposure_file <- exposure_files[i]
  cat(paste0("正在处理暴露文件 ", i, "/", length(exposure_files), ": ", exposure_file, "\n"))
  
  exp_dat <- read.csv(exposure_file)
  if (length(exp_dat) == 0) {
    cat("该暴露文件为空，跳过。\n")
    next  # 跳过本次循环，进入下一次循环
  }
  gene_name <- sub("^.*_", "", sub("\\.csv$", "", basename(exposure_file)))
  
  # 读取所有 GWAS 文件，提取所有 id
  all_ids <- c()
  cat("正在从所有 GWAS 文件中提取 ID...\n")
  for (gwas_file in gwas_files) {
    gwas_data <- read.table(gwas_file, header = FALSE, sep = "\t")
    all_ids <- c(all_ids, unique(gwas_data$V1))
  }
  all_ids <- unique(all_ids)
  cat(paste0("共提取到 ", length(all_ids), " 个唯一 ID。\n"))
  
  # 定义批量大小
  batch_size <- 10  # 可以根据实际情况调整
  num_batches <- ceiling(length(all_ids) / batch_size)
  
  for (batch_idx in 1:num_batches) {
    cat(paste0("正在处理第 ", batch_idx, "/", num_batches, " 个批次的 ID...\n"))
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(all_ids))
    batch_ids <- all_ids[start_idx:end_idx]
    
    # 提前构造结果文件路径并检查文件是否存在，根据 gene_name 和 id 的组合判断有效 id
    valid_ids <- c()
    for (id in batch_ids) {
      # 这里可以添加更复杂的判断逻辑，例如根据 gene_name 和 id 的组合是否满足某种条件
      # 目前简单假设只要 gene_name 和 id 组合有意义就认为是有效 id
      valid_ids <- c(valid_ids, id)
    }
    
    if (length(valid_ids) == 0) {
      cat("该批次没有有效 ID，跳过。\n")
      next  # 如果该批次没有有效 ID，跳过本次循环
    }
    cat(paste0("该批次有 ", length(valid_ids), " 个有效 ID 需要处理。\n"))
    
    # 并发请求 outcome_dat
    all_outcome_dat <- NULL
    cat("开始并发请求结果数据...\n")
    if (is_windows) {
      # 在 Windows 系统上使用 parLapply
      cl <- makeCluster(detectCores() - 1)
      # 设置集群节点输出打印信息到控制台
      clusterEvalQ(cl, options(warn = 1))
      clusterExport(cl, c("extract_outcome_data", "exp_dat"))
      all_outcome_dat <- parLapply(cl, valid_ids, function(id) {
        outcome_dat <- NULL
        
        # 初始化标志变量
        error_to_next <- FALSE
        while (TRUE) {
          result <- tryCatch({
            outcome_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = id)
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
          
          if (error_to_next == TRUE) {
            break
          }
          if (is.null(result)) {
            outcome_dat <- NULL
            break
          }
        }
        return(outcome_dat)
      })
      stopCluster(cl)
    } else {
      # 在 Unix 系统上使用 mclapply
      all_outcome_dat <- mclapply(valid_ids, function(id) {
        outcome_dat <- NULL
        
        # 初始化标志变量
        error_to_next <- FALSE
        while (TRUE) {
          result <- tryCatch({
            outcome_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = id)
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
          
          if (error_to_next == TRUE) {
            break
          }
          if (is.null(result)) {
            outcome_dat <- NULL
            break
          }
        }
        return(outcome_dat)
      }, mc.cores = detectCores() - 1)
    }
    cat("结果数据请求完成。\n")
    
    # 对每个有效 id 进行处理
    for (j in seq_along(valid_ids)) {
      id <- valid_ids[j]
      outcome_dat <- all_outcome_dat[[j]]
      
      # 构造结果文件路径
      dat_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_dat.csv"))
      res_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_res.csv"))
      or_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_or.csv"))
      
      # 数据匹配
      if (!is.null(outcome_dat)) {
        cat(paste0("正在处理有效 ID ", j, "/", length(valid_ids), ": ", id, "...\n"))
        dat <- harmonise_data(exp_dat, outcome_dat)
        
        if (!length(dat) || dim(dat)[1] == 0) {
          cat(paste0("ID ", id, " 数据匹配失败，跳过。\n"))
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
          cat(paste0("ID ", id, " 结果数据为空，跳过。\n"))
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
        cat(paste0("ID ", id, " 处理完成，结果已保存。\n"))
      }
    }
  }
  cat(paste0("暴露文件 ", exposure_file, " 处理完成。\n"))
}