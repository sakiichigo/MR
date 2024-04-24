remotes::install_github("WSpiller/RadialMR")

library(RadialMR)
#radialMR
egger_radial<- egger_radial(r_input = dat, alpha = 0.05,
                            weights = 1, summary = TRUE)

if(dim(dat)[1]>5){
  ivw_radial<- ivw_radial(r_input = dat, alpha = 0.05,
                          weights = 1, tol = 0.0001, summary = TRUE)
  #which(ivw_radial$data$Outliers == "Outlier") 
}
#filter
Outliers=ivw_radial$data[ivw_radial$data$Outliers=="Outlier",]$SNP

if(FALSE){
  #write
  write(Outliers,"C:/Users/user/Desktop/Outliers.txt")
  
  #read
  Outliers=data.table::fread("C:/Users/user/Desktop/Outliers.txt",header = FALSE)
  Outliers=Outliers$V1
}
#pick
exposure_dat=exposure_dat[!exposure_dat$SNP%in%Outliers,]