## ------------------------------------------------------------------------
library(tpAUC)

## ---- fig.height=5 , fig.width=5, fig.align='center'---------------------
library('pROC')
data(aSAH)
tproc.est(aSAH$outcome, aSAH$s100b, threshold=c(0.8,0.2),plot=T)
#estimate two-way partial AUC and generate a plot

## ----fig.height=5 , fig.width=5, fig.align='center'----------------------
proc.est(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8,plot=T)
# use method 'expect' to estimate partial AUC and draw a plot
proc.ci(aSAH$outcome,aSAH$s100b, cp=0.95 ,threshold=0.8,method='expect')
# use method 'expect' to infer partial AUC

## ------------------------------------------------------------------------
proc(aSAH$outcome,aSAH$s100b,threshold=0.8, method='expect',ci=TRUE, cp=0.95, plot=F)
# set ci=TRUE to get confidence interval

## ----fig.height=5 , fig.width=5, fig.align='center'----------------------
podc.est(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8,plot=T)
# estimate FNR partial ODC with method 'expect'
podc.ci(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8, cp=0.97)
# infer FNR partial ODC with method 'expect'

## ------------------------------------------------------------------------
podc(aSAH$outcome, aSAH$s100b,threshold=0.8, method='expect',ci=TRUE, cp=0.97, plot=F)
# inference and estimation

