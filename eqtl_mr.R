library(TwoSampleMR)
library(xlsx)
library(ieugwasr)

chunk_file_path <- "D:/document/bioInfo/alps-eye/allIdPath.txt"
if(TRUE){
  #F-Statistics
  Ff<-function(data){
    maf<-ifelse(data$eaf.exposure < 0.5, 
                data$eaf.exposure,
                1-data$eaf.exposure)
    n<-data$samplesize.exposure
    eaf<-data$eaf.exposure
    beta<-data$beta.exposure
    se<-data$se.exposure
    rr<-(2*(beta^2)*eaf*(1-eaf)) /(2*(beta^2)*eaf* (1-eaf) +2*n*eaf*(1-eaf)*se^2)
    if(all(is.na(rr))){
      rr<-beta^2/(beta^2+se^2*(data$samplesize.exposure-2))
    }
    ff<-((n-2)*rr)/(1-rr)
    mean<-mean(ff)
    li<-list(r2=mean(rr),
             #fs=ff,
             fm=mean)
    return(li)
  }
  
  #mr_raps_modified
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),silent = T)
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
  #eaf add
  snp_add_eaf <- function(dat, build = "37", pop = "EUR")
  {
    stopifnot(build %in% c("37","38"))
    stopifnot("SNP" %in% names(dat))
    
    # Create and get a url
    server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
    pop <- paste0("1000GENOMES:phase_3:",pop)
    
    snp_reverse_base <- function(x)
    {
      x <- stringr::str_to_upper(x)
      stopifnot(x %in% c("A","T","C","G"))
      switch(x,"A"="T","T"="A","C"="G","G"="C")
    }
    
    res_tab <- lapply(1:nrow(dat), function(i)
    {
      tryCatch({
        
        print(paste0("开始搜索第", i, "个SNP"))
        dat_i <- dat[i,]
        
        ext <- paste0("/variation/Homo_sapiens/",dat_i$SNP, "?content-type=application/json;pops=1")
        url <- paste(server, ext, sep = "")
        res <- httr::GET(url)
        
        # Converts http errors to R errors or warnings
        httr::stop_for_status(res)
        
        # Convert R objects from JSON
        res <- httr::content(res)
        res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
        
        # for no eaf
        if(length(res_pop)>0){
          # Filter query results based on population set
          res_pop <- try(res_pop[res_pop$population == pop,])
          if("try-error" %in% class(res_pop))
          {
            print(paste0("没有", pop, "群体的信息"))
            queried_effect_allele <- "NR"
            queried_other_allele <- "NR"
            queried_eaf <- -1
          }
          else
          {
            if(nrow(res_pop)==0)
            {
              print(paste0("没有", pop, "群体的信息"))
              queried_effect_allele <- "NR"
              queried_other_allele <- "NR"
              queried_eaf <- -1
            }
            else
            {
              queried_effect_allele <- res_pop[1,"allele"][[1]]
              queried_other_allele <- res_pop[2,"allele"][[1]]
              queried_eaf <- res_pop[1,"frequency"][[1]]    
            }
          }
          
          effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                                  dat_i$effect_allele.exposure,
                                  dat_i$effect_allele)
          
          other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                                 dat_i$other_allele.exposure,
                                 dat_i$other_allele)
          
          if("effect_allele.exposure" %in% names(dat))
          {
            name_output <- unique(c(names(dat), "eaf.exposure","reliability.exposure"))
          }
          else
          {
            name_output <- unique(c(names(dat), "eaf","reliability.exposure"))
          }
          
          len_effect_allele <- nchar(effect_allele)
          len_other_allele <- nchar(other_allele)
          
          if(len_effect_allele==1&len_other_allele==1)
          {
            if((queried_effect_allele==effect_allele & queried_other_allele==other_allele)|
               (queried_effect_allele==other_allele & queried_other_allele==effect_allele))
            {
              dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                           queried_eaf,
                                           1-queried_eaf)
              dat_i$eaf <- dat_i$eaf.exposure 
              dat_i$reliability.exposure <- "high"
            }
            else
            {
              r_queried_effect_allele <- snp_reverse_base(queried_effect_allele)
              r_queried_other_allele <- snp_reverse_base(queried_other_allele)
              if((r_queried_effect_allele==effect_allele & r_queried_other_allele==other_allele)|
                 (r_queried_effect_allele==other_allele & r_queried_other_allele==effect_allele))
              {
                dat_i$eaf.exposure <- ifelse(effect_allele == r_queried_effect_allele,
                                             queried_eaf,
                                             1-queried_eaf)
                dat_i$eaf <- dat_i$eaf.exposure 
                dat_i$reliability.exposure <- "high"
              }
              else
              {
                dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                             queried_eaf,
                                             1-queried_eaf)
                dat_i$eaf <- dat_i$eaf.exposure 
                dat_i$reliability.exposure <- "low"
              }
            }
          }
          
          else
          {
            # To identify the potential DEL/ INS
            short_allele <- ifelse(len_effect_allele==1,
                                   effect_allele,
                                   other_allele)
            short_allele_eaf <- ifelse(short_allele == queried_effect_allele, 
                                       queried_eaf, 
                                       1-queried_eaf)
            dat_i$eaf.exposure <- ifelse(effect_allele == short_allele,
                                         short_allele_eaf,
                                         1-short_allele_eaf)
            dat_i$eaf <- dat_i$eaf.exposure 
            dat_i$reliability.exposure <- "low"
          }
          
          dat_i[name_output]
          
        }
      }
      ,error=function(e) {
        print(e)
      })
      
    })
    
    return(do.call(rbind, res_tab))
  }
}
# 初始化错误数据向量
errorData_out <- c()
errorData_exp <- c()

# 结果保存目录
result_dir <- "D:/document/bioInfo/alps-eye/eqtl/result"
# 创建结果目录（如果不存在）
if (!dir.exists(result_dir)) {
  cat("结果目录不存在，正在创建：", result_dir, "\n")
  dir.create(result_dir, recursive = TRUE)
}

# 读取 D:/document/bioInfo/alps-eye/eqtl 目录下的 chunk_1.txt 文件
cat("正在读取文件：", chunk_file_path, "\n")
exposure_files <- readLines(chunk_file_path)

# 获取 GWAS 汇总数据目录下的所有 txt 文件
gwas_dir <- "D:/document/bioInfo/alps-eye/outcome"
cat("正在获取GWAS汇总数据目录下的txt文件：", gwas_dir, "\n")
gwas_files <- list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE)

outcomeId=c()
outPath <- list.files(path = gwas_dir, pattern = "\\.txt$", all.files = FALSE, full.names = FALSE)
cat("正在处理GWAS文件以获取outcomeId...\n")
for(m in outPath){
  readId=read.table(paste(gwas_dir,'/',m,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
  outcomeId=c(outcomeId,readId)
}
outcomeId=outcomeId[grepl("-", outcomeId)]
outcomeId=unique(outcomeId)
length_outcome=length(outcomeId)
cat("共获取到", length_outcome, "个唯一的outcomeId\n")

# 对读取到的文件进行循环
for (exposure_file in exposure_files) {
  cat("正在处理暴露文件：", exposure_file, "\n")
  if (!file.exists(exposure_file)) {
    cat("暴露文件不存在，跳过。\n")
    next  # 跳过本次循环，进入下一次循环
  }
  exp_dat <- read.csv(exposure_file)
  if (length(exp_dat) == 0) {
    cat("该暴露文件为空，跳过。\n")
    next  # 跳过本次循环，进入下一次循环
  }
  gene_name <- sub("^.*_", "", sub("\\.csv$", "", basename(exposure_file)))
  
  for (id in outcomeId) {
    # 构造结果文件路径
    dat_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_dat.csv"))
    res_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_res.csv"))
    or_file <- file.path(result_dir, paste0("eqtl_", gene_name, "_", id, "_or.csv"))
    
    # 检查文件是否存在，如果存在则跳过本次循环
    if (file.exists(dat_file) && file.exists(res_file) && file.exists(or_file)) {
      cat("结果文件已存在，跳过：", dat_file, res_file, or_file, "\n")
      next
    }
    
    # 初始化标志变量
    message_to_next <- TRUE
    error_to_next <- FALSE
    
    cat("正在尝试获取outcome数据，outcomeId为：", id, "\n")
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
          message("检测到达到OpenGWAS配额限制，等待5分钟后继续...")
          Sys.sleep(300) # 等待5分钟（300秒）
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
      cat("开始进行数据匹配...\n")
      dat <- harmonise_data(exp_dat, outcome_dat)
      
      if (!length(dat) || dim(dat)[1] == 0) {
        cat("数据匹配失败，跳过。outcomeId: ", id, " SNP: ", exp_dat$SNP, "\n")
        errorData_out <- c(errorData_out, id)
        errorData_exp <- c(errorData_exp, exp_dat$SNP)
        next
      }
      
      # 异质性多效性
      cat("正在计算异质性多效性...\n")
      ple <- mr_pleiotropy_test(dat)
      he <- mr_heterogeneity(dat)
      
      if (length(ple$pval) == 0 || length(he$Q_pval) == 0) {
        # 可能出现没有足够SNP
        cat("可能没有足够SNP，正在进行MR分析...\n")
        res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
        ivw <- "Inverse variance weighted"
      } else {
        if (he[he$method == "Inverse variance weighted", ]$Q_pval < 0.05) {
          cat("存在异质性，更改方法进行异质性分析和MR分析...\n")
          he <- mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre")) # he更改方法
          # 存在异质性时随机效应ivw_mre，同时raps替换simple mode
          res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_raps_modified", "mr_weighted_mode"))
          ivw <- "Inverse variance weighted (multiplicative random effects)"
        } else {
          cat("不存在异质性，进行MR分析...\n")
          res <- mr(dat, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_raps_modified", "mr_weighted_mode"))
          ivw <- "Inverse variance weighted"
        }
      }
      
      # 空数据
      if (nrow(res) == 0) {
        cat("结果数据为空，跳过。\n")
        next
      }
      
      res[is.na(res$method), ]$method <- "Robust Adjusted Profile Score"
      or <- generate_odds_ratios(res)
      
      # or写入HE，只保留IVW
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
      
      # or写入PLE
      if (length(ple) == 0) {
        or[, 16] <- NA
      } else {
        or[, 16] <- ple[1, 7]
      }
      colnames(or)[16] <- "ple"
      
      or[, 17] <- Ff(exp_dat)$fm
      colnames(or)[17] <- "F"
      
      cat("正在保存结果文件...\n")
      write.csv(dat, dat_file, row.names = FALSE)
      write.csv(res, res_file, row.names = FALSE)
      write.csv(or, or_file, row.names = FALSE)
      cat("结果文件保存成功：", dat_file, res_file, or_file, "\n")
    }
  }
}