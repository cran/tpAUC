% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/podc.est.R
\name{podc.est}
\alias{podc.est}
\title{Partial ODC Estimation}
\usage{
podc.est(response, predictor, threshold = 0.9, method = "MW",
  smooth = FALSE)
}
\arguments{
\item{response}{a factor, numeric or character vector of responses; 
typically encoded with 0 (negative) and 1 (positive). 
Only two classes can be used in a ROC curve. If its levels are not 0 and 1,
the first level will be defaultly regarded as negative.}

\item{predictor}{a numeric vector of the same length than response, containing the predicted value of each observation. An ordered factor is coerced to a numeric.}

\item{threshold}{numeric; false negative rate (FNR) constraint.}

\item{method}{methods to estimate partial ODC. \code{MW}: Mann-Whitney statistic. \code{expect}: method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} adapted from \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.}

\item{smooth}{if \code{TRUE}, the ODC curve is passed to \code{\link[pROC]{smooth}} to be smoothed.}
}
\value{
Estimation of FNR partial ODC.
}
\description{
Estimate the area of region under ordinal dominance curve with pre-specific FNR constraint (FNR-pODC). See \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} for details.
}
\details{
This function estimates FNR partial ODC given response, predictor and pre-specific FNR constraint.  
         \code{MW}: Mann-Whitney statistic. \code{expect}: method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017} adapted from \href{http://www.ncbi.nlm.nih.gov/pubmed/20729218}{Wang and Chang, 2011}. \code{jackknife}: jackknife method in \href{http://www3.stat.sinica.edu.tw/statistica/j27n1/j27n117/j27n117.html}{Yang et al., 2017}.
}
\examples{

library('pROC')
data(aSAH)
podc.est(aSAH$outcome, aSAH$s100b, method='expect',threshold=0.8 )

}
\author{
Hanfang Yang, Kun Lu, Xiang Lyu, Feifang Hu, Yichuan Zhao.
}
\seealso{
\code{\link[tpAUC]{proc.est}}
}

