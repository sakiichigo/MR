#install.packages("LDlinkR")
library(LDlinkR) 
library(xlsx)

#https://ldlink.nih.gov/?tab=ldlink

#exposure
length=length(exposure_dat$SNP)
LDinfo=data.frame()
num=floor(length/10)
if(length>10){
  for(i in 1:num){
    snp=exposure_dat[(i*10-9):(i*10),]$SNP
    LDinfo_ <- LDtrait(snps = snp, 
                       pop = "EUR", r2d = "r2", 
                       r2d_threshold = 0.001,
                       win_size=10000,
                       token = '0e18a58e5aea', 
                       file ="FALSE")
    LDinfo=rbind(LDinfo,LDinfo_)
  }
}
snp=exposure_dat[((num+1)*10-9):length,]$SNP
LDinfo_ <- LDtrait(snps = snp, 
                   pop = "EUR", r2d = "r2", 
                   r2d_threshold = 0.001,
                   win_size=10000,
                   token = '0e18a58e5aea', 
                   file ="FALSE")
LDinfo=rbind(LDinfo,LDinfo_)
#filter A|B|C
snp_filter=unique(LDinfo[grep("diabetic r|ocular pressure|laucoma|ataract",LDinfo$GWAS_Trait),]$Query)

if(FALSE){
  #write
  name=unique(exposure_dat$id.exposure)
  write(snp_filter,paste0("C:/Users/user/Desktop/",name," snp_filter.txt"))
  write.xlsx(LDinfo,paste0("C:/Users/user/Desktop/",name," LDinfo.xlsx"))
  
  #read
  snp_filter=data.table::fread("C:/Users/user/Desktop/snp_filter.txt",header = FALSE)
  snp_filter=snp_filter$V1
}
#pick
exposure_dat=exposure_dat[exposure_dat$SNP%in%snp_filter,]

