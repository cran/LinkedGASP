\name{GASP_plot}
\alias{GASP_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the GASP
}
\description{
Function allows to plot the GASP in case of one-dimensional input.
}
\usage{
GASP_plot(em, fun, data, emul_type, labels, yax, ylab, xlab,ylim,
col_CI_area,col_points,col_fun,col_mean,plot_training = FALSE, plot_fun = TRUE)
}
\arguments{
  \item{em}{
the returned output from the function eval_type1_GASP(...) or eval_type2_GASP(...).
}
  \item{fun}{
Simulator function. Currently only one-dimensional input is supported.
}
  \item{data}{
Training data and smoothness. The same as supplied to eval_GASP_RFP(...) for construction of the GASP.
}
  \item{emul_type}{
A text string which provides description of an emulator.
}
  \item{labels}{
As in standard R plot.
}
  \item{yax}{
As in standard R plot.
}
  \item{ylab}{
As in standard R plot.
}
  \item{xlab}{
As in standard R plot.
}
  \item{ylim}{
As in standard R plot.
}
  \item{col_CI_area}{
Color of a credible area.
}
  \item{col_points}{
Color of the training points.
}
  \item{col_fun}{
Color of a simulator function.
}
  \item{col_mean}{
Color of the emulator of the GASP mean.
}
  \item{plot_training}{
(Not) to plot the training points. Default is FALSE. 
}
\item{plot_fun}{
(Not) to plot the simulator function. Default is TRUE. 
}
}
\value{
Plot
}
\author{
Ksenia N. Kyzyurova, kseniak.ucoz.net
}
\note{
The function requires further development to be automated for visualization along a single dimension out of multiple dimensions and along two dimensions out of multiple dimensions.
}
\examples{
## Function f1 is a simulator
f1<-function(x){sin(pi*x)}

## One-dimensional inputs are x1
x1 <- seq(-1,1,.37)

## The following contains the list of data inputs (training) and outputs (fD) together with the 
## assumed fixed smoothness of a computer model output.
data.f1 <- list(training = x1,fD = f1(x1), smooth = 1.99)

## Evaluation of GASP parameters
f1_MLEs = eval_GASP_RFP(data.f1,list(function(x){x^0},function(x){x^1}),1,FALSE)

## Evaluate the emulator
xn = seq(-1,1,.01)
GASP_type1_f1 <- eval_type1_GASP(as.matrix(xn),f1_MLEs)

## Plot the emulator
par(mfrow = c(1,1))
par(mar = c(6.1, 6.1, 5.1, 2.1))
ylim = c(-1.5,1.5)
GASP_plot(GASP_type1_f1,fun = f1,data = data.f1,"",ylim = ylim, plot_training = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }% use one of  RShowDoc("KEYWORDS")
