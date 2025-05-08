# 加载所需的包
library(TwoSampleMR)
library(VariantAnnotation)

# 定义暴露数据所在目录
exposure_dir <- "E:/eqtl"
# 获取 E:/eqtl 目录下的所有 csv 文件
exposure_files <- list.files(exposure_dir, pattern = "\\.csv$", full.names = TRUE)

# 定义 GWAS 数据所在目录
gwas_dir <- "E:\\sacropeniaGwas"
# 获取 E:\sacropeniaGwas 目录下的所有 .vcf.gz 文件
gwas_files <- list.files(gwas_dir, pattern = "\\.vcf.gz$", full.names = TRUE)

# 定义结果保存目录
result_dir <- "D:/document/bioInfo/heart-sarco-0330/eqtl/result2"

# 初始化错误数据向量
errorData_out <- c()
errorData_exp <- c()

# 对 E:/eqtl 目录下的 csv 文件进行循环
for (exposure_file in exposure_files) {
  
  exp_dat <- read.csv(exposure_file)
  
  # 如果 exp_dat 为空，跳过本次循环
  if (length(exp_dat) == 0) {
    next
  }
  
  # 循环处理每个 .vcf.gz 文件
  for (vcfFileName in gwas_files) {
    result <- tryCatch({
      # 读取 VCF 文件
      vcf <- readVcf(vcfFileName, genome = "hg19")  # 假设基因组版本为 hg19，可根据实际情况修改
      
      # 提取所需信息
      snps <- rowRanges(vcf)$ID
      # 使用 ref 和 alt 函数提取等位基因信息
      effect_allele <- as.character(ref(vcf))
      other_allele <- sapply(alt(vcf), function(x) as.character(x)[1])
      
      # 这里需要根据实际 VCF 文件中的信息来计算或提取 beta、se、pval、eaf 等信息
      # 以下是示例，假设这些信息存储在 INFO 字段中，需要根据实际情况修改
      beta <- info(vcf)$A1_beta
      se <- info(vcf)$se
      pval <- info(vcf)$pval
      eaf <- info(vcf)$A1_freq
      
      # 创建 outcome_dat 数据框
      outcome_dat <- data.frame(
        SNP = snps,
        A1 = effect_allele,
        A2 = other_allele,
        A1_beta = beta,
        se = se,
        pval = pval,
        A1_freq = eaf
      )
      
      # 筛选出与 exp_dat 中 SNP 匹配的行
      outcome_dat <- outcome_dat[outcome_dat$SNP %in% exp_dat$SNP, ]
    },
    error = function(e) {
      print(e)
      return(NULL)
    })
    
    if (is.null(result)) {
      next
    }
    
    if (!is.null(outcome_dat)) {
      dat <- harmonise_data(exp_dat, outcome_dat)
      
      if (!length(dat) || dim(dat)[1] == 0) {
        errorData_out <- c(errorData_out, vcfFileName)
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
      
      # 保存结果
      id <- tools::file_path_sans_ext(basename(vcfFileName))
      dat_file <- file.path(result_dir, paste0("eqtl_", id, "_dat.csv"))
      res_file <- file.path(result_dir, paste0("eqtl_", id, "_res.csv"))
      or_file <- file.path(result_dir, paste0("eqtl_", id, "_or.csv"))
      
      write.csv(dat, dat_file, row.names = FALSE)
      write.csv(res, res_file, row.names = FALSE)
      write.csv(or, or_file, row.names = FALSE)
    }
  }
}