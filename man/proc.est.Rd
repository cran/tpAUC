% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proc.est.R
\name{proc.est}
\alias{proc.est}
\title{Partial AUC Estimation}
\usage{
proc.est(response, predictor, threshold = 0.9, method = "MW",
  smooth = FALSE)
}
\arguments{
\item{response}{a factor, numeric or character vector of responses; 
typically encoded with 0 (negative) and 1 (positive). 
Only two classes can be used in a ROC curve. If its levels are not 0 and 1,
the first level will be defaultly regarded as negative.}

\item{predictor}{a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.}

\item{threshold}{numeric; false positive rate (FPR) constraint.}

\item{method}{methods to estimate FPR-pAUC. \code{MW}: Mann-Whitney statistic. \code{expect}: method in (2.2) \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.}

\item{smooth}{if \code{TRUE}, the ROC curve is passed to \code{\link[pROC]{smooth}} to be smoothed.}
}
\value{
Estimate of FPR partial AUC.
}
\description{
Estimate the area of region under ROC curve with pre-specific FPR constraint (FPR-pAUC). See \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} for details.
}
\details{
This function estimates FPR partial AUC given response, predictor and pre-specific FPR constraint.
         \code{MW}: Mann-Whitney statistic. \code{expect}: method in (2.2) \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
}
\examples{

library('pROC')
data(aSAH)
proc.est(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8)

}
\author{
Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
}
\seealso{
\code{\link[tpAUC]{tproc.est}}, \code{\link[tpAUC]{podc.est}}
}

