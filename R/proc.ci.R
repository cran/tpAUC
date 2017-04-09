#' Partial AUC Inference
#'
#' Infer the area of region under ROC curve with pre-specific FPR constraint (FPR-pAUC). See \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} for details.
#' 
#' @param response a factor, numeric or character vector of responses; 
#'    typically encoded with 0 (negative) and 1 (positive). 
#'    Only two classes can be used in a ROC curve. If its levels are not 0/1,
#'    the first level will be defaultly regarded as negative.
#'    
#' @param predictor a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.
#' @param cp numeric; coverage probability of confidence interval. 
#' @param threshold numeric; false positive rate (FPR) constraint.
#' @param method methods to estimate FPR-pAUC. \code{MW}: Mann-Whitney statistic. \code{expect}: method in (2.2) \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' 
#' @details This function infers FPR partial AUC given response, predictor and pre-specific FPR constraint. 
#'          \code{MW}: Mann-Whitney statistic. method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} adapted from \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' @return Confidence interval of FPR partial AUC. 
#'
#' @author Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
#' @seealso \code{\link[tpAUC]{tproc.est}}, \code{\link[tpAUC]{podc.ci}}
#'
#' @examples
#' 
#' library('pROC')
#' data(aSAH)
#' proc.ci(aSAH$outcome, aSAH$s100b, cp=0.95 ,threshold=0.9,method='expect')
#' 
#' @export 
#'
#' @import pROC
#' @importFrom  stats ecdf qnorm
#'          
#'        
proc.ci=function(response,predictor, cp=0.95 ,threshold=0.9,method='MW') {
  
  if ( any(is.na(response) | is.na(predictor))){
    warning('NA will be removed.')
    # remove NA
    nonna=(!is.na(response)) & (!is.na(predictor))
    response=response[nonna]
    predictor=predictor[nonna]
  } 
  
  if  (any(levels(response)!=c(0,1)) |( length(levels(response))!=2) ) {
    warning('response levels are not 0/1, the first level is defaultly regarded as negative.')
  }
  
  if (!is.numeric(predictor)) {
    warning('predictor is coerced to a numeric.')
    predictor=as.numeric(predictor)
  }
  
  
  if (length(response) != length(predictor)) {
    stop('response and predictor should have the same length.')
  } else if ( any( (threshold > 1) | (threshold < 0) | (length(threshold)!=1)) ){ 
    stop('threshold should between 0 and 1.')
  } else if (!any( method == c('MW', 'expect', 'jackknife')))  {
    stop('please input a correct method type.')
  } else if ( any((length(cp)!=1) | (cp > 1) | (cp < 0)) ) {
    stop('cp should be a numeric between 0 and 1.')
  }
  
  
  
  level=levels(as.factor(response))
  # the levels of response
  neg_level=level[1];pos_level=level[2]
  # default the first level is negative 
  
  neg_predictor=predictor[response==neg_level]
  pos_predictor=predictor[response==pos_level]
  # extract negative and positice predictor
  
  FPR=threshold
  # thresholding FPR in the ROC graph
  
  neg_size=length(neg_predictor)
  pos_size=length(pos_predictor)
  # sample size 
  
  #quantile of negaitve/positive sample
  
  
  es= 1- ecdf(neg_predictor)(pos_predictor) # empirical survival function
  es[es>FPR]=FPR
  pauctilde=FPR - sum(es)/pos_size    
  
  pauch=c()
  
  for (h in 1:pos_size) {
    es= (1- ecdf(neg_predictor)(pos_predictor))[-h] # empirical survival function
    es[es>FPR]=FPR
    pauch[h]=FPR - mean(es)
  }
  
  for ( h in (pos_size+1):(pos_size+neg_size)) {
    Sgnh=c()
    for( i in 1:pos_size) { 
      Sgnh[i]= min( mean(neg_predictor[-(h-pos_size)] > pos_predictor[i]) , FPR)
    }
    
    pauch[h]=FPR -  mean(Sgnh)
    
  }
  
  Vh=(pos_size+neg_size)*pauctilde-(pos_size+neg_size-1)*pauch
  mu=mean(Vh)
  S2= mean((Vh- mu)^2)
  
  
  if (method=='MW') {
    mu=proc.est(response,predictor,threshold=threshold,method='MW' ,smooth=FALSE)
    ci=c(mu-qnorm(1- (1-cp)/2) * sqrt(S2)/sqrt(pos_size+neg_size) , mu- qnorm((1-cp)/2) *sqrt(S2)/sqrt(pos_size+neg_size))
  } else if (method=='expect') {
    mu=proc.est(response,predictor,threshold=threshold,method='expect' ,smooth=FALSE)
    ci=c(mu-qnorm(1- (1-cp)/2) * sqrt(S2)/sqrt(pos_size+neg_size) , mu- qnorm((1-cp)/2) *sqrt(S2)/sqrt(pos_size+neg_size))
  } else if (method=='jackknife') { 
    ci=c(mu-qnorm(1- (1-cp)/2) * sqrt(S2)/sqrt(pos_size+neg_size) , mu- qnorm((1-cp)/2) *sqrt(S2)/sqrt(pos_size+neg_size))
  }
  ci=matrix(ci,byrow=T,1,2)
  colnames(ci)=c(paste(((1-cp)/2)*100, '%'),paste((1- (1-cp)/2)*100, '%'))
  return(ci)
}