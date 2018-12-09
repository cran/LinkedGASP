\name{TGASPmetrics}
\alias{TGASPmetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performance measurement of a T-GASP
}
\description{
Evaluates frequentist performance of a T-GASP.
}
\usage{
TGASPmetrics(TGASP, true_output, ref_output)
}
\arguments{
  \item{TGASP}{
TGASP emulator (in the paper this is done within an objective Bayesian implementation - OB emulator.)
}
  \item{true_output}{
  Output from the simulator.
}
  \item{ref_output}{
Heuristic emulator output.
}
}
\details{
See examples which illustrate the use of the function.
}
\value{
List of performance measures.
\item{RMSPE_base }{Root mean square predictive error with respect to the heuristic emulator output.}
\item{RMSPE }{Root mean square predictive error for the emulator output}
\item{ratio }{ratio of RMSPE_base to RMSPE. Ratio = RMSPE_base/RMSPE}
\item{CIs }{95\% central credible intervals}
\item{emp_cov }{95\% empirical coverage within the CIs}
\item{length_CIs }{Average lenght of 95\% central credible intervals}
}
\references{Ksenia N. Kyzyurova, James O. Berger, and Robert L. Wolpert. Coupling computer models through linking their statistical emulators. SIAM/ASA Journal on Uncertainty Quantification, 6(3): 1151-1171, 2018}
\author{
Ksenia N. Kyzyurova, ksenia.ucoz.net
}
\examples{
## Function f1 is a simulator
f1<-function(x){sin(pi*x)}

## One-dimensional inputs are x1
x1 <- seq(-1,1,.37)

## The following contains the list of data inputs (training) and outputs (fD) together with 
## the assumed fixed smoothness of a computer model output.
data.f1 <- list(training = x1,fD = f1(x1), smooth = 1.99)

## Evaluation of GASP parameters
f1_MLEs = eval_GASP_RFP(data.f1,list(function(x){x^0},function(x){x^1}),1,FALSE)

## Evaluate the emulator
xn = seq(-1,1,.01)
TGASP_f1 <- eval_TGASP(as.matrix(xn),f1_MLEs)

## Plot the emulator
par(mfrow = c(1,1))
par(mar = c(6.1, 6.1, 5.1, 2.1))
ylim = c(-1.5,1.5)
TGASP_plot(TGASP_f1,f1,data.f1,ylim = ylim)

## Measure the performance of the emulator
TGASPmetrics(TGASP_f1,f1(xn),mean(f1(xn)))
}
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
