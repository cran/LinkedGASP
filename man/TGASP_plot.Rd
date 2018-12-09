\name{TGASP_plot}
\alias{TGASP_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
T-GASP plot
}
\description{
Function allows to plot the TGASP in case of one-dimensional input. Black-and-white version.}
\usage{
TGASP_plot(tem, fun, data, labels, ylim, points)
}
\arguments{
  \item{tem}{
TGasP emulator.
}
  \item{fun}{
Simulator function.
}
\item{data}{
Training data and smoothness. The same as supplied to eval_GASP_RFP(...) for construction of a GASP.
}
  \item{labels}{
As in standard R plot.
}
  \item{ylim}{
As in standard R plot.
}
  \item{points}{
(Not) to plot the training points.
}
}
\details{
See examples.
}
\value{
Plot
}
\author{
Ksenia N. Kyzyurova, kseniak.ucoz.net
}
\note{
The function requires further development to be automated for visualization along a single dimension out of multiple dimensions and along two dimensions out of multiple dimensions.

This function needs to be automated to allow for fast visualization of a single emualtor (with no comparison to the actual simulator function), etc.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }% use one of  RShowDoc("KEYWORDS")
