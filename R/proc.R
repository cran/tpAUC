#' Partial AUC Estimation and Inference
#'
#' Estimate and infer the area of region under ROC curve with pre-specific FPR constraint (FPR-pAUC). See \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} for details.
#' 
#' @param response a factor, numeric or character vector of responses; 
#'    typically encoded with 0 (negative) and 1 (positive). 
#'    Only two classes can be used in a ROC curve. If its levels are not 0/1,
#'    the first level will be defaultly regarded as negative.
#'    
#' @param predictor a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.
#' @param threshold numeric; false positive rate (FPR) constraint.
#' @param method methods to estimate FPR-pAUC. \code{MW}: Mann-Whitney statistic. \code{expect}: method in (2.2) \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' @param ci logic; compute the confidence interval of estimation?
#' @param cp numeric; coverage probability of confidence interval. 
#' @param smooth if \code{TRUE}, the ROC curve is passed to \code{\link[pROC]{smooth}} to be smoothed.
#' 
#' @details This function estimates and infers FPR partial AUC given response, predictor and pre-specific FPR constraint.
#'          \code{MW}: Mann-Whitney statistic. \code{expect}: method in (2.2) \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
#' @return Estimate and Inference of FPR partial AUC. 
#'
#' @author Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
#' @seealso \code{\link[pROC]{roc}}, \code{\link[tpAUC]{tproc.est}}, \code{\link[tpAUC]{proc.est}}, \code{\link[tpAUC]{proc.ci}}
#'
#' @examples
#' 
#' library('pROC')
#' data(aSAH)
#' proc(aSAH$outcome, aSAH$s100b,threshold=0.9, method='expect',ci=TRUE, cp=0.95)
#' 
#' @export 
#'
#' @import pROC
#' @importFrom   stats  ecdf approxfun runif
#' @importFrom graphics lines abline      
#'          
#'        
proc=function(response,predictor,threshold=0.9, method='MW',ci=TRUE, cp=0.95, smooth=FALSE) {
  pauc=proc.est(response,predictor,threshold=threshold, method=method,smooth=smooth)
  l=list(pauc=pauc)
  if (ci==TRUE){
    c=proc.ci(response,predictor, cp=cp ,threshold=threshold,method=method )
    l[['ci']]=c
  }
  return(l)
}
  