library(TwoSampleMR)
library(xlsx)
library(ggplot2)
library(ieugwasr)
library(MRPRESSO)
library(RadialMR)
#Define folders
workPath="C:/Users/user/Desktop/work"

#Do you keep or he ple records
saveAllRecord=TRUE
exportPDF=TRUE

#Create folders
dir.create(paste(workPath,"/outcome",sep=""))#outcome file path
dir.create(paste(workPath,"/exposure",sep=""))#exposure file path
#Put the exposure and ending IDs in /exposure and /outcome
dir.create(paste(workPath,"/result",sep=""))#result file path
if(saveAllRecord){
  dir.create(paste(workPath,"/result/or",sep=""))#or file path
  dir.create(paste(workPath,"/result/ple",sep=""))#ple file path
  dir.create(paste(workPath,"/result/he",sep=""))#he file path
  dir.create(paste(workPath,"/result/pdf",sep=""))#pdf file path
  dir.create(paste(workPath,"/result/steiger",sep=""))#steiger file path
  dir.create(paste(workPath,"/result/direct",sep=""))#direct file path
  dir.create(paste(workPath,"/result/presso",sep=""))#presso file path
  dir.create(paste(workPath,"/result/radialMR",sep=""))#radialMR file path
}

#The following one-click runs

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
  #rr<-beta^2/(beta^2+se^2*(data$samplesize.exposure-2))
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


#outcome file
outcomeId=c()
outPath=list.files(path=paste(workPath,"/outcome",sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
for(m in outPath){
  readId=read.table(paste(workPath,'/outcome/',m,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
  outcomeId=c(outcomeId,readId)
}
outcomeId=unique(outcomeId)
length_outcome=length(outcomeId)
#Expose files + analyze + save
expPath=list.files(path=paste(workPath,"/exposure",sep=""), pattern=NULL, all.files=FALSE, full.names=FALSE)
length_exporsure_path=length(expPath)
errorData_exp=c();#error dataa
errorData_out=c();
k_temp=1#operation progress
i_temp=1
j_temp=1
#An error is interrupted during the long execution cycle, and the K cycle can be continued if the cycle is re-executed (the environment is not cleared).
for(k in k_temp:length_exporsure_path){
  k_temp=k
  sheet_name=expPath[k]
  exposureId=read.table(paste(workPath,'/exposure/',sheet_name,sep=""),fill = TRUE,row.names = NULL)[,c(1)]
  length_exporsure=length(exposureId)
  resultTable=data.frame()
  resultTable_ivw=data.frame()
  for(i in i_temp:length_exporsure){#exposure
    i_temp=i
    n=exposureId[i]
    while(TRUE){
      message_to_next <<- TRUE
      error_to_next <<- FALSE
      try({
        withCallingHandlers(
          exposure_dat<-extract_instruments(n),
          message = function(c) if (stringr::str_detect(as.character(c),"Failed to"))
            message_to_next <<- FALSE)
        error_to_next <<- TRUE})
      if(message_to_next == TRUE&error_to_next == TRUE) { break }
    }
    #empty data
    if(!length(exposure_dat) ){
      errorData_exp=c(errorData_exp,n)
      next
    }
    #outcome
    for(j in j_temp:length_outcome){
      j_temp=j
      m=outcomeId[j]
      while(TRUE){
        message_to_next <<- TRUE
        error_to_next <<- FALSE
        try({
          withCallingHandlers(
            outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,outcomes = m),
            message = function(c) if (stringr::str_detect(as.character(c),"Failed to"))
              message_to_next <<- FALSE)
          error_to_next <<- TRUE})
        if(message_to_next == TRUE&error_to_next == TRUE) { break }
      }
      #empty data
      if(!length(outcome_dat) ){
        errorData_out=c(errorData_out,m)
        next
      }
      dat <- harmonise_data(exposure_dat,outcome_dat)
      he<-mr_heterogeneity(dat);
      if(he[he$method=="Inverse variance weighted",]$Q_pval<0.05){
        he=mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw_mre"))#he change method
        #Random effects when heterogeneity is present: ivw_mre
        res <- mr(dat,method_list=c( "mr_egger_regression","mr_weighted_median","mr_ivw_mre","mr_raps_modified","mr_weighted_mode"))
        ivw="Inverse variance weighted (multiplicative random effects)"
      }else{
        #res <- mr(dat)
        res <- mr(dat,method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_raps_modified","mr_weighted_mode"))
        ivw="Inverse variance weighted"
      }
      res[is.na(res$method),]$method="Robust Adjusted Profile Score"
      #mr_method_list()#List the methods
      #empty data
      if(length(res)==0){
        next
      }
      or<-generate_odds_ratios(res);
      ple<-mr_pleiotropy_test(dat);
      print(as.POSIXlt(Sys.time()))
      print(paste('exposure file:',sheet_name,'   process:',k,'/',length_exporsure_path))
      print(paste('id.exposure:',n,'   process:',i,'/',length_exporsure))
      print(paste('id.outcome:',m,'   process:',j,'/',length_outcome))
      #Record keeping
      if(saveAllRecord){
        write.csv(or,paste(workPath,"/result/or/",n,' ',m," or.csv",sep=""),row.names = F)
        write.csv(ple,paste(workPath,"/result/ple/",n,' ',m," ple.csv",sep=""),row.names = F)
        write.csv(he,paste(workPath,"/result/he/",n,' ',m," he.csv",sep=""),row.names = F)
        if((!is.na(dat$samplesize.exposure[1]))&(!is.na(dat$samplesize.outcome[1]))){
          direct=directionality_test(dat)#Reverse cause and effect
          steiger<-steiger_filtering(dat)#Steiger filter
          write.csv(direct,paste(workPath,"/result/direct/",n,' ',m," direct.csv",sep=""),row.names = F)
          write.csv(steiger,paste(workPath,"/result/steiger/",n,' ',m," steiger.csv",sep=""),row.names = F)
        }
        presso=mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                         SignifThreshold = 0.05)
        write.csv(presso$`Main MR results`,paste(workPath,"/result/presso/",n,' ',m," Main MR results.csv",sep=""),row.names = F)
        
        write.csv(presso$`MR-PRESSO results`,paste(workPath,"/result/presso/",n,' ',m," MR-PRESSO results.csv",sep=""),row.names = F)
        #radialMR
        egger_radial<- egger_radial(r_input = dat, alpha = 0.05,
                                    weights = 1, summary = TRUE)
        
        write.csv(egger_radial$data,paste(workPath,"/result/radialMR/",n,' ',m," egger_radial data.csv",sep=""),row.names =FALSE)
        if(dim(dat)[1]>5){
          ivw_radial<- ivw_radial(r_input = dat, alpha = 0.05,
                                  weights = 1, tol = 0.0001, summary = TRUE)
          write.csv(ivw_radial$data,paste(workPath,"/result/radialMR/",n,' ',m," ivw_radial data.csv",sep=""),row.names =FALSE)
          #which(ivw_radial$data$Outliers == "Outlier") 
          plot_radial(c(ivw_radial,egger_radial))
        }else{
          plot_radial(egger_radial)
        }
        ggsave(paste0(workPath,"/result/pdf/",n,' ',m," radialMR.pdf"))
        
      }
      if(exportPDF){
        #leave out
        single <- mr_leaveoneout(dat)
        if(length(unique(single$SNP))==length(single$SNP)){
          mr_leaveoneout_plot(single)
          ggsave(paste0(workPath,"/result/pdf/",n,' ',m," leaveoneout.pdf"))
        }
        #scatter
        mr_scatter_plot(res,dat)
        ggsave(paste0(workPath,"/result/pdf/",n,' ',m," scatter_plot.pdf"))
        res_single <- mr_singlesnp(dat)
        if(length(unique(single$SNP))==length(single$SNP)){
          #forest
          mr_forest_plot(res_single)
          ggsave(paste0(workPath,"/result/pdf/",n,' ',m," forest_plot.pdf"))
          #funnel
          mr_funnel_plot(res_single)
          ggsave(paste0(workPath,"/result/pdf/",n,' ',m," funnel_plot.pdf"))
        }
      }
      #select results (save if one is significant)
      pval=or$pval
      if(length(pval)>1){
        for(p in pval)
          if(p<0.05){
            #or writes to he, leaving only ivw
            if(dim(he)[1]>1&he[2,]$method==ivw){
              or[,15]=he[he$method==ivw,]$Q_pval
            }else{
              or[,15]=he[1,]$Q_pval
            }
            colnames(or)[15]="he"
            #or write to ple
            or[,16]=ple[1,7]
            colnames(or)[16]="ple"
            or[,17]=Ff(exposure_dat)$fm
            colnames(or)[17]="F"
            #resultTable: Summarize five methods, one of which is significant
            resultTable=rbind(resultTable,or)
            #resultTable_ivw: Further sensitivity checks and IVW methods only
            if(dim(or[or$method==ivw,])[1]>0&or[or$method==ivw,]$pval<0.05&length(ple$pval)>0){
              if(ple[1,]$pval>0.05){
                if(length(he$method)==2&he[2,]$method==ivw&he[2,]$Q_pval>0.05){
                  resultTable_ivw=rbind(resultTable_ivw,or[or$method==ivw,])
                }else if(length(he$method)==1&he[1,]$method==ivw&he[1,]$Q_pval>0.05){
                  resultTable_ivw=rbind(resultTable_ivw,or[or$method==ivw,])
                }
              }
            }
            break;
          }
      }
    }
    j_temp=1
  }
  i_temp=1
  #save
  if(length(resultTable)!=0){
    write.xlsx(resultTable,paste(workPath,'/result/result.xlsx',sep=""), sheetName=strsplit(sheet_name,"\\.",fixed = FALSE)[[1]][1],append=TRUE,row.names=FALSE) 
    if(length(resultTable_ivw)!=0){
      write.xlsx(resultTable_ivw,paste(workPath,'/result/result-ivw.xlsx',sep=""), sheetName=strsplit(sheet_name,"\\.",fixed = FALSE)[[1]][1],append=TRUE,row.names=FALSE) 
    }
    print(paste('file:',sheet_name,' saved '))
  }else{
    print(paste('file:',sheet_name,' length 0'))
  }
}
print('errorData exposure:')
unique(errorData_exp)
print('errorData outcome:')
unique(errorData_out)



