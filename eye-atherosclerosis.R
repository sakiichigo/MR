library(TwoSampleMR)
library(xlsx)
library(ggplot2)
library(ieugwasr)
library(MRPRESSO)
library(RadialMR)
##options(ieugwasr_api = 'https://gwas-api.mrcieu.ac.uk/')
options(ieugwasr_api = 'https://api.opengwas.io/api/')
###1.读取数据
#方法
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
        
        print(paste0("seaching for No.", i, " SNP"))
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
            print(paste0("There is not information for population ",pop))
            queried_effect_allele <- "NR"
            queried_other_allele <- "NR"
            queried_eaf <- -1
          }
          else
          {
            if(nrow(res_pop)==0)
            {
              print(paste0("There is not information for population ",pop))
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
#用exp作暴露
workPath="D:/document/bioInfo/atherosclerosis-eye"
dir.create(workPath)


outPath="E:/anaFiles/atherosclerosisGwas/tsv/"

#导出or,he,ple记录
saveRecord=TRUE
#导出presso,radialMR,plot文件 (开启循环变慢)
exportFile=FALSE

#Part2:需要添加文件部分
#创建文件夹
dir.create(paste(workPath,"/outcome",sep=""))#结局文件路径
dir.create(paste(workPath,"/exposure",sep=""))#暴露文件路径
#把暴露和结局ID放进/exposure和/outcome
dir.create(paste(workPath,"/result",sep=""))#结果路径
if(saveRecord){
  dir.create(paste(workPath,"/result/or",sep=""))#or文件路径
  dir.create(paste(workPath,"/result/ple",sep=""))#ple文件路径
  dir.create(paste(workPath,"/result/he",sep=""))#he文件路径
  dir.create(paste(workPath,"/result/pdf",sep=""))#pdf文件路径
  dir.create(paste(workPath,"/result/steiger",sep=""))#steiger文件路径
  dir.create(paste(workPath,"/result/direct",sep=""))#direct文件路径
  dir.create(paste(workPath,"/result/presso",sep=""))#presso文件路径
  dir.create(paste(workPath,"/result/radialMR",sep=""))#radialMR文件路径
}

errorData_exp=c();#报错数据
errorData_out=c();
allResultTable=data.frame()
resultTable=data.frame()
resultTable_ivw=data.frame()
#结局文件


outLists=list.files(path=outPath, pattern=NULL, all.files=FALSE, full.names=FALSE)
#
outfilename=outLists[3]
#expLists=expLists[which(regexpr(thisName,expLists)==1):length(expLists)]
#开始分析
outcomeId=c()
outcomePath="/outcome"
outPath_=list.files(path=paste(workPath,outcomePath,sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
for(m in outPath_){
  readId=read.table(paste(workPath,outcomePath,'/',m,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
  outcomeId=c(outcomeId,readId)
}
outcomeId=outcomeId[grepl("-", outcomeId)]
outcomeId=unique(outcomeId)
length_outcome=length(outcomeId)


for(n in outcomeId){
  while(TRUE){
    message_to_next <<- TRUE
    error_to_next <<- FALSE
    result=tryCatch({
      withCallingHandlers(
        exposure_dat<-extract_instruments(n),
        message = function(c) if (stringr::str_detect(as.character(c),"Failed to"))
          message_to_next <<- FALSE)
      error_to_next <<- TRUE
    },
    error=function(e) {
      print(e)
      if (grepl("used up your OpenGWAS allowance", as.character(e), fixed = TRUE)) {
        message("检测到达到 OpenGWAS 配额限制，等待 5 分钟后继续...")
        Sys.sleep(300) # 等待 5 分钟（300 秒）
      }
      return(NULL)
    })
    if(is.null(result)){
      break
    }
    if(message_to_next == TRUE&error_to_next == TRUE)
      break
  }
  if(is.null(result)){
    next
  }
  #空数据
  if(!length(exposure_dat)||dim(exposure_dat)[1]==0 ){
    errorData_exp=c(errorData_exp,n)
    next
  }
  thisName=strsplit(outfilename,"\\.",fixed = FALSE)[[1]][1]
  m=thisName
  outCsvName=paste0(outPath,thisName,".h.tsv")
  result=tryCatch({
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = outCsvName,
      sep = "\t",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "standard_error",
      eaf_col = "effect_allele_frequency",
      effect_allele_col = "effect_allele", 
      other_allele_col = "other_allele",
      pval_col = "p_value",
      pos="base_pair_location"
      #eaf_col = "EFFECTALLELEFREQUENCY"
    )
  },
  error=function(e) {
    print(e)
    return(NULL)
  })
  if(is.null(result)){
    next
  }
  
  
  #空数据
  if(!length(outcome_dat)||dim(outcome_dat)[1]==0 ){
    print(paste0(m," is null"))
    errorData_out=c(errorData_out,m)
    next
  }
  outcome_dat$id.outcome=thisName
  dat <- harmonise_data(exposure_dat,outcome_dat)
  
  if(!length(dat)||dim(dat)[1]==0 ){
    errorData_out=c(errorData_out,m)
    errorData_exp=c(errorData_exp,n)
    next
  }
  #异质性多效性
  ple<-mr_pleiotropy_test(dat);
  he<-mr_heterogeneity(dat);
  if(length(ple$pval)==0||length(he$Q_pval)==0){#可能出现没有足够snp
    
    #res <- mr(dat)
    res <- mr(dat,method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_raps_modified","mr_weighted_mode"))
    if(nrow(res)==0)
      next
    ivw="Inverse variance weighted"
    if(!is.null(res$method))
      res[is.na(res$method),]$method="Robust Adjusted Profile Score"
  }else{
    if(he[he$method=="Inverse variance weighted",]$Q_pval<0.05){
      he=mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre"))#he更改方法
      #存在异质性时随机效应ivw_mre，同时raps替换simple mode
      res <- mr(dat,method_list=c( "mr_egger_regression","mr_weighted_median","mr_ivw_mre","mr_raps_modified","mr_weighted_mode"))
      if(nrow(res)==0)
        next
      ivw="Inverse variance weighted (multiplicative random effects)"
    }else{
      #res <- mr(dat)
      res <- mr(dat,method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_raps_modified","mr_weighted_mode"))
      if(nrow(res)==0)
        next
      ivw="Inverse variance weighted"
    }
    if(!is.null(res$method))
      res[is.na(res$method),]$method="Robust Adjusted Profile Score"
  }
  #mr_method_list()#列出方法
  #空数据
  if(length(res)==0){
    next
  }
  or<-generate_odds_ratios(res);
  #or写入HE，只保留IVW
  if(length(he)==0){
    or[,15]=NA
  }else{
    if(dim(he)[1]>1&he[2,]$method==ivw){
      or[,15]=he[he$method==ivw,]$Q_pval
    }else{
      or[,15]=he[1,]$Q_pval
    }
  }
  colnames(or)[15]="he"
  #or写入PLE
  if(length(ple)==0){
    or[,16]=NA
  }else{
    or[,16]=ple[1,7]
  }
  colnames(or)[16]="ple"
  
  or[,17]=Ff(exposure_dat)$fm
  colnames(or)[17]="F"
  allResultTable=rbind(allResultTable,or)
  
  print(as.POSIXlt(Sys.time()))
  print(paste('id.exposure(i):',n))
  print(paste('id.outcome(j):',m))
  
  #留存记录
  if(saveRecord){
    write.csv(or,paste(workPath,"/result/or/",n,' ',m," or.csv",sep=""),row.names = F)
    write.csv(ple,paste(workPath,"/result/ple/",n,' ',m," ple.csv",sep=""),row.names = F)
    write.csv(he,paste(workPath,"/result/he/",n,' ',m," he.csv",sep=""),row.names = F)
    if((!is.null(dat$samplesize.exposure[1]))&(!is.null(dat$samplesize.outcome[1]))){
      if((!is.na(dat$samplesize.exposure[1]))&(!is.na(dat$samplesize.outcome[1]))){
        direct=directionality_test(dat)#反向因果
        steiger<-steiger_filtering(dat)#Steiger过滤
        write.csv(direct,paste(workPath,"/result/direct/",n,' ',m," direct.csv",sep=""),row.names = F)
        write.csv(steiger,paste(workPath,"/result/steiger/",n,' ',m," steiger.csv",sep=""),row.names = F)
      }
    }
    
    
  }
  
  
}
