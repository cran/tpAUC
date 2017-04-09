#' Partial ODC Estimation
#'
#' Estimate the area of region under ordinal dominance curve with pre-specific FNR constraint (FNR-pODC). See \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} for details.
#' 
#' @param response a factor, numeric or character vector of responses; 
#'    typically encoded with 0 (negative) and 1 (positive). 
#'    Only two classes can be used in a ROC curve. If its levels are not 0 and 1,
#'    the first level will be defaultly regarded as negative.
#' @param predictor a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.
#' @param threshold numeric; false negative rate (FNR) constraint.
#' @param method methods to estimate partial ODC. \code{MW}: Mann-Whitney statistic. \code{expect}: method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} adapted from \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' @param smooth if \code{TRUE}, the ODC curve is passed to \code{\link[pROC]{smooth}} to be smoothed.
#' 
#' @details This function estimates FNR partial ODC given response, predictor and pre-specific FNR constraint.  
#'          \code{MW}: Mann-Whitney statistic. \code{expect}: method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} adapted from \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' @return Estimation of FNR partial ODC. 
#'
#' @author Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
#' @seealso \code{\link[tpAUC]{proc.est}}
#'
#' @examples
#' 
#' library('pROC')
#' data(aSAH)
#' podc.est(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8 )
#' 
#' @export 
#'
#' @import pROC
#' @importFrom    stats ecdf approxfun runif
#' @importFrom graphics lines abline
#'        
podc.est=function(response,predictor,threshold=0.9,method='MW' ,smooth=FALSE) {
  
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
  }
  
  
  
  level=levels(as.factor(response))
  # the levels of response
  neg_level=level[1];pos_level=level[2]
  # default the first level is negative 
  
  neg_predictor=predictor[response==neg_level]
  pos_predictor=predictor[response==pos_level]
  # extract negative and positice predictor
  
  FNR=threshold
  # thresholding FNR in the ROC graph
  
  neg_size=length(neg_predictor)
  pos_size=length(pos_predictor)
  # sample size 
  
  #quantile of negaitve/positive sample
  
  
  if (method=='MW') {
    
    if (floor(FNR*pos_size)==0) {
      neg_quant=sort(pos_predictor)[1];
    } else {
      neg_quant=sort(pos_predictor)[floor(FNR*pos_size)]  
    }
    
    
    V=c();l=0
    # pairwise indicator
    
    for (i in 1:pos_size) {
      for (j in 1:neg_size) {
        l=l+1
        V[l]= ( ( !( pos_predictor[i] < neg_predictor[j])) &  ( ! ( pos_predictor[i] > neg_quant ) ) )
        # pairwise comparison
      }
    }
    
    podc = mean(V)  # estimated partial odc
    
  } else if (method=="expect") {
    
    es= ecdf(pos_predictor)(neg_predictor) # empirical distribution function
    es[es>FNR]=FNR
    podc=FNR - mean(es)
    
  } else if (method=='jackknife') {
    
    es= ecdf(pos_predictor)(neg_predictor) # empirical distribution function
    es[es>FNR]=FNR
    podctilde = FNR - mean(es)
    
    podch=c()
    
    for (h in 1:neg_size) {
      es= (ecdf(pos_predictor)(neg_predictor))[-h] # empirical survival function
      es[es>FNR]=FNR
      podch[h]=FNR - mean(es)
    }
    
    for ( h in (neg_size+1):(pos_size+neg_size)) {
      Gnh=c()
      for( i in 1:neg_size) { 
        Gnh[i]= min( mean(  !(pos_predictor[-(h-neg_size)] >  neg_predictor[i]) ) , FNR)
      }
      
      podch[h]=FNR -  mean(Gnh)
      
    }
    
    Uh=(pos_size+neg_size)*podctilde-(pos_size+neg_size-1)*podch
    podc=mean(Uh)
  }
  
 
  return(podc)
  
}

