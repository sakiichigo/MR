library(MendelianRandomization)
library(MVMR)
library(TwoSampleMR)
library(ieugwasr)
id_exposure <- c("ukb-b-12440","finn-b-DM_RETINOPATHY","ukb-b-12397","finn-b-H7_GLAUCOMA_XFG")

id_outcome <- "ukb-a-379"

while(TRUE){
  message_to_next <<- TRUE
  error_to_next <<- FALSE
  try({
    withCallingHandlers(
      exposure_dat <- mv_extract_exposures(id_exposure) ,
      message = function(c) if (stringr::str_detect(as.character(c),"Failed to"))
        message_to_next <<- FALSE)
    error_to_next <<- TRUE})
  if(message_to_next == TRUE&error_to_next == TRUE) { break }
}
while(TRUE){
  message_to_next <<- TRUE
  error_to_next <<- FALSE
  try({
    withCallingHandlers(
      outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) ,
      message = function(c) if (stringr::str_detect(as.character(c),"Failed to"))
        message_to_next <<- FALSE)
    error_to_next <<- TRUE})
  if(message_to_next == TRUE&error_to_next == TRUE) { break }
}

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat) 
res <- mv_multiple(mvdat) 
res 


SummaryStats<-cbind(mvdat[["outcome_beta"]],
                    mvdat[["exposure_beta"]][,1],
                    mvdat[["exposure_beta"]][,2],
                    mvdat[["exposure_beta"]][,3],
                    mvdat[["exposure_beta"]][,4],
                    mvdat[["exposure_se"]][,1],
                    mvdat[["exposure_se"]][,2],
                    mvdat[["exposure_se"]][,3],
                    mvdat[["exposure_se"]][,4],
                    mvdat[["outcome_se"]])
SummaryStats<-data.frame(SummaryStats)

presso=mr_presso(BetaOutcome ="X1", BetaExposure = c("X2","X3","X4","X5"), SdOutcome ="X10", 
                 SdExposure = c("X6","X7","X8","X9"), 
                 OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  
                 SignifThreshold = 0.05)

MRMVInputObject_1<-mr_mvinput(
  bx=cbind(SummaryStats$X2,SummaryStats$X3,SummaryStats$X4,SummaryStats$X5),
  bxse=cbind(SummaryStats$X6,SummaryStats$X7,SummaryStats$X8,SummaryStats$X9),
  by = SummaryStats$X1, 
  byse = SummaryStats$X10)
MRMVInputObject_1

mvivw_=data.frame()
mvegger_=data.frame()
mvmedian_=data.frame()
  
#ivw
MRMVObject_ivw <- mr_mvivw(MRMVInputObject_1, 
                           model = "default",
                           correl = FALSE,
                           distribution = "normal",
                           alpha = 0.05)
MRMVObject_ivw

mvivw_=rbind(mvivw_,data.frame(MRMVObject_ivw$Exposure,MRMVObject_ivw$Estimate,MRMVObject_ivw$StdError,
                                                 MRMVObject_ivw$CILower,MRMVObject_ivw$CIUpper,MRMVObject_ivw$Pvalue))
#egger
MRMVObject_egger<-mr_mvegger(MRMVInputObject_1,
                             orientate = 1,
                             correl = FALSE,
                             distribution = "normal",
                             alpha = 0.05)
MRMVObject_egger

mvegger_=rbind(mvegger_,data.frame(MRMVObject_egger$Exposure,MRMVObject_egger$Estimate,MRMVObject_egger$StdError.Est,
                                   MRMVObject_egger$CILower.Est,MRMVObject_egger$CIUpper.Est,MRMVObject_egger$Pvalue.Est))
#median
MRMVObject_median<-mr_mvmedian(
  MRMVInputObject_1,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject_median

mvmedian_=rbind(mvmedian_,data.frame(MRMVObject_median$Exposure,MRMVObject_median$Estimate,MRMVObject_median$StdError,
                                     MRMVObject_median$CILower,MRMVObject_median$CIUpper,MRMVObject_median$Pvalue))
###MVMR
F.data <- format_mvmr(BXGs = SummaryStats[,c(2,3,4,5)],
                      BYG = SummaryStats[,1],
                      seBXGs = SummaryStats[,c(6,7,8,9)],
                      seBYG = SummaryStats[,10],
                      RSID = rownames(SummaryStats))

head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
ivwMvmrRes <- ivw_mvmr(r_input = F.data)

write.xlsx(mvivw_,"C:/Users/user/Desktop/mvivw_.xlsx")
write.xlsx(mvegger_,"C:/Users/user/Desktop/mvegger_.xlsx")
write.xlsx(mvmedian_,"C:/Users/user/Desktop/mvmedian_.xlsx")
