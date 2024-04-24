#install.packages("meta")
library("meta")

#read

lnor<- log(or[,"or"])
lnuci<- log(or[,"or_uci95"])
lnlci<- log(or[,"or_lci95"])
selnor<- (lnuci-lnlci)/(2*1.96)#95%CI

or=data.table::fread("C:/Users/user/Desktop/cataract.csv")
or=as.data.frame(or)
for (i in 1:length(or$exposure)){
  or[i,]$exposure=strsplit(or[i,]$exposure,"\\|",fixed = FALSE)[[1]][1]
}
#or=or[or$id.exposure!="",]
#Meta
pfs=metagen(studlab=or$id.exposure,TE=or$be,seTE  =or$se,lower = or$lo_ci,upper = or$up_ci,
            #overall = F,#common model
            #random = F,#random model
            sm="or")

pfs
title="dr"
write.csv(pfs,paste(workPath,"/result/meta/",title," meta.csv",sep=""),row.names =FALSE)



pfs$TE.fixed
pfs$lower.fixed
pfs$upper.fixed
pfs$pval.fixed

#forest


pdf(file=paste0("C:/Users/user/Desktop/hgs_final/result/forest/",title," forest_plot.pdf"),width=14,height=28)
forest(pfs);
dev.off()


