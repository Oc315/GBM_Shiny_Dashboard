# This code is translated from Jing Qiu's Matlab code ---xq 11/1/05

# This code is to do equivalence analysis  for two sample comparison 

# test H_0: |\theta_g|^2/max(\sigma_g^2, \sigma_00^2) \geq \Delta
#      H_a: |\theta_g|^2/max(\sigma_g^2, \sigma_00^2) < \Delta


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Canonical model:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#D_g~ N(\theta_g, c\sigam_g^2)
#MSE_g ~ \sigma^2_g \chi^2_g/df 
#sigma_D=\sqrt(C)\sigma_g
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The statistical prcedure(estimation method) is
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# if \sigma_g^2<\sigma_00^2; |D_g|<\sqrt{\Delta}\sigma_00-t\sigma_D(TOST)
# if \sigma_g^2>=\sigma_00^2; |D_g|<\sigma_D \sqrt{F}  (UMPI)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# and the corresponding p-value is 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# p1=P(F<(|D|^2/(c\sigma_g^2)) (pvalue for UMPI test
# p2=P(T>= (\sqrt(Delta\sigma_00-|D|)/(\sigma_g\sqrt(c))) (p-value for TOST)
# p=p2 if \sigma_g^2<\sigma_00^2
#  =p1 if \sigma_g^2>\sigma_00^2 
#=min(p1,p2) if \sigma_g^2=\sigma_00^2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%the inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# filename is a textfile with three columns: sample mean for group 1;
#        sample mean for group 2 and the residual variance for the gene
# c   is the scale that connects the variance of the estimator and the
#  residuals; sigma_D^2=csigma_g^2, where sigma_g^2 is the residual
#  variance and sigma^2_D is the varaince of the difference in VG.
# df is the degrees of freedom for the residual variance
# sigma00 is the constant in the hypothesis
# Delta is the threshold for the hypothsis.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
equivtest2comp <- function (filename,c, df, sigma00, Delta)

{
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in the raw data: A(:1) is the sample mean for group 1 and 
# A(:,2) is the sample mean for group 2.A(:,3) is the residaul variance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read in the VG and residual variance estimator
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
A=filename
  
D=A[,1]-A[,2] #D_g~ N(\theta_g, \sigam_g^2(1/n1+1/n2))
var0=A[,3] # residual variance
sigma0=sqrt(var0)
varD=var0*c#the variance for the difference
sigmaD=sqrt(varD)

#delta=sqrt(Delta)*sigma00
delta=log2(Delta)

ncpF=(2*Delta^2-1)/c

#%%%%%%%%%%%calculate the p-value for UMPI test
obs1=(abs(D)/sigmaD)^2
pval1=pf(obs1,df1=1,df2=df,ncp=ncpF, lower.tail=T)
indx1=which(sigma0>sigma00)

#%%%%%%%%%calculate the p-value for the TOST test
obs2=(delta-abs(D))/sigmaD
pval2=1-pt(obs2,df, lower.tail=T) 
indx2=which(sigma0<sigma00)

#3%%%%%%%%%calculate the p-value of our proposed test
indx3=which(sigma0==sigma00)
pval=pval1
pval[indx2]=pval2[indx2]
pval3=apply(cbind(pval1,pval2), 1, min)
if (is.null(indx3)==0){
  pval[indx3]=pval3[indx3]
}
out=cbind(pval1, pval2, pval)
out
}
#the output is the p-value for all the genes against three


