#' Two-Way Partial AUC Estimation
#'
#' Estimate the area of region under ROC curve under pre-specific FPR/TPR constraints (two-way partial AUC). See \href{http://arxiv.org/abs/1508.00298}{Yang et al., 2016} for details.
#' 
#' @param response a factor, numeric or character vector of responses; 
#'    typically encoded with 0 (negative) and 1 (positive). 
#'    Only two classes can be used in a ROC curve. If its levels are not 0 and 1,
#'    the first level will be defaultly regarded as negative.
#'    
#' @param predictor a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.
#' @param threshold a length-two numeric vector; the first element is FPR threshold, the second is TPR.  
#' @param smooth if \code{TRUE}, the ROC curve is passed to \code{\link[pROC]{smooth}} to be smoothed.
#' 
#' @details This function estimates two-way partial AUC given response, predictor and pre-specific FPR/TPR constraints. 
#'
#' @return Estimate of two-way partial AUC. 
#'
#' @author Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
#' @seealso \code{\link[pROC]{roc}}, \code{\link[tpAUC]{podc.est}}, \code{\link[tpAUC]{proc.est}}
#'
#' @examples
#' 
#' library('pROC')
#' data(aSAH)
#' tproc.est(aSAH$outcome, aSAH$s100b, threshold=c(0.8,0.2))
#' 
#'
#' @export 
#'
#' @import pROC
#' @importFrom   stats  approxfun runif
#' @importFrom graphics lines abline
#'          
tproc.est=function(response, predictor, threshold=c(1,0), smooth=FALSE) {
  
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
  } else if ( any( (threshold > 1) | (threshold < 0) ) ){ 
    stop('TPR/FPR should between 0 and 1.')
  } 
  
    
  level=levels(as.factor(response))
  # the levels of response
  neg_level=level[1];pos_level=level[2]
  # default the first level is negative 
  
  neg_predictor=predictor[response==neg_level]
  pos_predictor=predictor[response==pos_level]
  # extract negative and positice predictor
  
  FPR=threshold[1];TPR=threshold[2]
  # thresholding in the ROC graph
  
  neg_size=length(neg_predictor)
  pos_size=length(pos_predictor)
  # sample size 
  
  #quantile of negaitve/positive sample
  
  if (floor((1-TPR)*pos_size) == 0 ) { # if FPR =1, quantile = 0 
    pos_quant=sort(pos_predictor)[1]
  } else {pos_quant=sort(pos_predictor)[floor((1-TPR)*pos_size)]}
  
  if (floor((1-FPR)*neg_size)==0) {
    neg_quant=sort(neg_predictor)[1];
  } else {
    neg_quant=sort(neg_predictor)[floor((1-FPR)*neg_size)]  
  }
  
  
  V=c();l=0
  # pairwise indicator
  
  for (i in 1:pos_size) {
    for (j in 1:neg_size) {
      l=l+1
      V[l]= ( ( !( pos_predictor[i] < neg_predictor[j])) & (! (pos_predictor[i] > pos_quant)) & ( ! (neg_predictor[j] < neg_quant ) ) )
      # pairwise comparison
    }
  }
  tpauc=sum(V)/(pos_size*neg_size) # estimated two-way partial AUC
 
  
  return(tpauc)
}


