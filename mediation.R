
results=data.frame()
#置信区间
library(RMediation)
beta_all=or$se#A to C(得到总效应，beta_all)
beta1=or1$b#A to B(得到beta1)
beta2=or2$b#B to C(得到beta2)
beta12=beta1*beta2
beta_dir=beta_all-beta12
se1=or1$se
se2=or2$se
medciRes=medci(mu.x=beta1,mu.y=beta2,se.x=se1,se.y=se2,rho=0,alpha=0.1,type="prodclin")

result=c(
  medciRes$Estimate,
  medciRes$`90% CI`,
  beta12,
  beta_dir,
  medciRes$SE
)

results=  rbind(results,result)
colnames(results)=c('Estimate','90% CI_Low','90% CI_Up','beta12','beta_dir','medciRes$SE')

write.xlsx(results,"C:/Users/user/Desktop/tsmr.xlsx")

