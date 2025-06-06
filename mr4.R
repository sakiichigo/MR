library(TwoSampleMR)
library(xlsx)
library(ggplot2)
library(ieugwasr)
library(MRPRESSO)
library(RadialMR)
#options(ieugwasr_api = 'https://gwas-api.mrcieu.ac.uk/')

#Part1:需要修改的部分
#定义文件夹
rootPath='D:/document/bioInfo'
workPath="D:/document/bioInfo/heart-sarco-0330"
dir.create(workPath)
#导出or,he,ple记录
getGwasInfo=FALSE
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

#Part3:一键执行部分

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


#结局文件
outcomeId=c()
if(FALSE){
  outcomePath="/outcome4"
  outPath=list.files(path=paste(workPath,outcomePath,sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
  for(m in outPath){
    readId=read.table(paste(workPath,outcomePath,'/',m,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
    outcomeId=c(outcomeId,readId)
  }
}

outcomeId=readId=read.table(paste(rootPath,'/liver/output/merged.txt',sep=""),fill = TRUE,row.names = NULL)[,c(1)]
outcomeId=outcomeId[grepl("-", outcomeId)]
outcomeId=unique(outcomeId)
length_outcome=length(outcomeId)

#暴露文件
expPath=list.files(path=paste(workPath,"/exposure",sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
length_exposure_path=length(expPath)
#开始分析
allResultTable=data.frame()
errorData_exp=c();#报错数据
errorData_out=c();
k_temp=1#运行进度，每层循环从外到内k>i>j (exposure_file>exposure_Id>outcome_Id)
i_temp=1
j_temp=1
#选择运行位置: 执行长循环时报错中断，将j_temp替换成j(或者跳过这一步j+1)，重新执行下面的循环可接续运行(k,i同理)
for(k in k_temp:length_exposure_path){#读取暴露id
  k_temp=k
  sheet_name=expPath[k]
  exposureId=read.table(paste(workPath,'/exposure/',sheet_name,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
  length_exposure=length(exposureId)
  resultTable=data.frame()
  resultTable_ivw=data.frame()
  for(i in i_temp:length_exposure){#读取暴露
    i_temp=i
    n=exposureId[i]
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
    #结局变量
    for(j in j_temp:length_outcome){
      print(as.POSIXlt(Sys.time()))
      print(paste('exposure file(k):',sheet_name,'   process:',k,'/',length_exposure_path))
      print(paste('id.exposure(i):',n,'   process:',i,'/',length_exposure))
      print(paste('id.outcome(j):',m,'   process:',j,'/',length_outcome))
      j_temp=j
      m=outcomeId[j]
      while(TRUE){
        message_to_next <<- TRUE
        error_to_next <<- FALSE
        result=tryCatch({
          withCallingHandlers(
            outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,outcomes = m),
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
      if(!length(outcome_dat)||dim(outcome_dat)[1]==0 ){
        errorData_out=c(errorData_out,m)
        next
      }
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
        ivw="Inverse variance weighted"
      }else{
        if(he[he$method=="Inverse variance weighted",]$Q_pval<0.05){
          he=mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre"))#he更改方法
          #存在异质性时随机效应ivw_mre，同时raps替换simple mode
          res <- mr(dat,method_list=c( "mr_egger_regression","mr_weighted_median","mr_ivw_mre","mr_raps_modified","mr_weighted_mode"))
          ivw="Inverse variance weighted (multiplicative random effects)"
        }else{
          #res <- mr(dat)
          res <- mr(dat,method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_raps_modified","mr_weighted_mode"))
          ivw="Inverse variance weighted"
        }
      }
      #mr_method_list()#列出方法
      #空数据
      if(nrow(res)==0){
        next
      }
      res[is.na(res$method),]$method="Robust Adjusted Profile Score"
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
      #写入gwas信息
      if(getGwasInfo){
        gwasinfo_exp=ieugwasr:: gwasinfo(n)
        gwasinfo_out=ieugwasr:: gwasinfo(m)
        if(is.null(gwasinfo_exp$sample_size)){
          if(!is.null(gwasinfo_exp$ncase)&!is.null(gwasinfo_exp$ncontrol)){
            gwasinfo_exp$sample_size=gwasinfo_exp$ncase+gwasinfo_exp$ncontrol
            exposure_dat$samplesize.exposure=gwasinfo_exp$sample_size
          }
        }
        if(is.null(gwasinfo_out$sample_size)){
          if(!is.null(gwasinfo_out$ncase)&!is.null(gwasinfo_out$ncontrol)){
            gwasinfo_out$sample_size=gwasinfo_out$ncase+gwasinfo_out$ncontrol
            outcome_dat$samplesize.outcome=gwasinfo_out$sample_size
          }
        }
        
        or[,18]=paste0(gwasinfo_exp$population,"|",gwasinfo_out$population)
        colnames(or)[18]="population"
        or[,19]=paste0(gwasinfo_exp$sample_size,"|",gwasinfo_out$sample_size)
        colnames(or)[19]="sample_size"
        or[,20]=paste0(gwasinfo_exp$nsnp,"|",gwasinfo_out$nsnp)
        colnames(or)[20]="nsnps"
        or[,21]=paste0(gwasinfo_exp$year,"|",gwasinfo_out$year)
        colnames(or)[21]="year"
        or[,22]=paste0(gwasinfo_exp$pmid,"|",gwasinfo_out$pmid)
        colnames(or)[22]="pmid"
        or[,23]=paste0(gwasinfo_exp$consortium,"|",gwasinfo_out$consortium)
        colnames(or)[23]="consortium"
      }
      allResultTable=rbind(allResultTable,or)
      
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
        if(exportFile){
          #可视化
          result=tryCatch({
            #留一
            single <- mr_leaveoneout(dat)
            if(length(unique(single$SNP))==length(single$SNP)){
              mr_leaveoneout_plot(single)
              ggsave(paste0(workPath,"/result/pdf/",n,' ',m," leaveoneout.pdf"), device = cairo_pdf, dpi = 300,width = 12, height = 8)
            }
            #散点
            mr_scatter_plot(res,dat)
            ggsave(paste0(workPath,"/result/pdf/",n,' ',m," scatter_plot.pdf"), device = cairo_pdf, dpi = 300,width = 12, height = 8)
            res_single <- mr_singlesnp(dat)
            if(length(unique(single$SNP))==length(single$SNP)){
              #森林
              mr_forest_plot(res_single)
              ggsave(paste0(workPath,"/result/pdf/",n,' ',m," forest_plot.pdf"), device = cairo_pdf, dpi = 300,width = 12, height = 8)
              #漏斗
              mr_funnel_plot(res_single)
              ggsave(paste0(workPath,"/result/pdf/",n,' ',m," funnel_plot.pdf"), device = cairo_pdf, dpi = 300,width = 12, height = 8)
            }
          },
          error=function(e) {
            print(e)
          })
          
          #敏感性校验
          #presso
          result=tryCatch({
            presso=mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                             SignifThreshold = 0.05)
            write.csv(presso$`Main MR results`,paste(workPath,"/result/presso/",n,' ',m," Main MR results.csv",sep=""),row.names = F)
            outliersIndices=presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
            if(length(outliersIndices)>=2){
              presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`=toString(outliersIndices)
            }
            write.csv(presso$`MR-PRESSO results`,paste(workPath,"/result/presso/",n,' ',m," MR-PRESSO results.csv",sep=""),row.names = F)
          },
          error=function(e) {
            print(e)
          })
          
          #radialMR
          result=tryCatch({
            egger_radial<- egger_radial(r_input = dat, alpha = 0.05,
                                        weights = 1, summary = TRUE)
            write.csv(egger_radial$data,paste(workPath,"/result/radialMR/",n,' ',m," egger_radial data.csv",sep=""),row.names =FALSE)
            ivw_radial<- ivw_radial(r_input = dat, alpha = 0.05,
                                    weights = 1, tol = 0.0001, summary = TRUE)
            write.csv(ivw_radial$data,paste(workPath,"/result/radialMR/",n,' ',m," ivw_radial data.csv",sep=""),row.names =FALSE)
            #which(ivw_radial$data$Outliers == "Outlier") 
            plot_radial(c(ivw_radial,egger_radial))
            ggsave(paste0(workPath,"/result/pdf/",n,' ',m," radialMR.pdf"), device = cairo_pdf, dpi = 300,width = 12, height = 8)
          },
          error=function(e) {
            print(e)
          })
          
        }
        
      }
      #挑选结果（其一显著则保存）
      pval=or$pval
      if(length(pval)>1){
        for(p in pval){
          if (is.na(p)) {
            next  # 如果是缺失值，跳过当前循环
          }
          if(p<0.05){
            #resultTable汇总五种方法其一显著
            resultTable=rbind(resultTable,or)
            #resultTable_ivw进一步进行敏感性校验且只含ivw方法
            # 假设 ivw 是一个已经定义的变量
            # 首先检查 or 数据框中 method 列等于 ivw 的行是否存在
            if (any(or$method == ivw)) {
              # 提取 method 列等于 ivw 的行
              ivw_rows <- or[or$method == ivw, ]
              # 检查 pval 列是否存在且没有缺失值
              if (!is.null(ivw_rows$pval) && all(!is.na(ivw_rows$pval))) {
                # 检查 ple 数据框的 pval 列是否存在且没有缺失值
                if (!is.null(ple$pval) && length(ple$pval) > 0 && all(!is.na(ple$pval))) {
                  if (dim(ivw_rows)[1] > 0 && ivw_rows$pval < 0.05) {
                    if (ple[1, ]$pval > 0.05) {
                      # 检查 he 数据框的情况
                      if (length(he$method) == 2 && he[2, ]$method == ivw && he[2, ]$Q_pval > 0.05) {
                        resultTable_ivw <- rbind(resultTable_ivw, ivw_rows)
                      } else if (length(he$method) == 1 && he[1, ]$method == ivw && he[1, ]$Q_pval > 0.05) {
                        resultTable_ivw <- rbind(resultTable_ivw, ivw_rows)
                      }
                    }
                  }
                }
              }
            }
            break;
          }
        }
          
      }
    }
    j_temp=1
  }
  i_temp=1
  #保存
#  if(length(resultTable)!=0){
 #   write.xlsx(resultTable,paste(workPath,'/result/result.xlsx',sep=""), sheetName=strsplit(sheet_name,"\\.",fixed = FALSE)[[1]][1],append=TRUE,row.names=FALSE) 
  #  if(length(resultTable_ivw)!=0){
   #   write.xlsx(resultTable_ivw,paste(workPath,'/result/result-ivw.xlsx',sep=""), sheetName=strsplit(sheet_name,"\\.",fixed = FALSE)[[1]][1],append=TRUE,row.names=FALSE) 
    #}
    #print(paste('file:',sheet_name,' saved '))
  #}else{
  #  print(paste('file:',sheet_name,' length 0'))
  #}
}
write.xlsx(allResultTable,paste(workPath,'/result/allResult.xlsx',sep=""), sheetName="all",append=TRUE,row.names=FALSE) 
if(FALSE){
  #遇到报错导致循环中断可执行此部分手动生成AllResult
  allResultTable=data.frame()
  getPath=list.files(path=paste(workPath,"/result/or",sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
  for(m in getPath){
    readCsv=read.csv(paste(workPath,'/result/or/',m,sep=""),fill = TRUE,row.names = NULL)
    allResultTable=rbind(allResultTable,readCsv)
  }
  write.xlsx(allResultTable,paste(workPath,'/result/allResult.xlsx',sep=""), sheetName="all",append=TRUE,row.names=FALSE) 
}
print('errorData exposure:')
unique(errorData_exp)
print('errorData outcome:')
unique(errorData_out)
#resultTable_ivw: 只含ivw方法并且无异质多效性汇总
#resultTable: 方法其一显著汇总
#allResultTable:  所有产生数据的汇总

