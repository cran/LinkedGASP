\name{eval_TGASP}
\alias{eval_TGASP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
T-GASP emulator
}
\description{
This function evaluates the third GASP of a computer model within objective Bayesian (OB) implementation of the GASP, resulting in T-GASP.}
\usage{
eval_TGASP(input, GASPparams)
}
\arguments{
  \item{input}{
Input values (the same dimension as training input data in the next argument GASPparams)
}
  \item{GASPparams}{
The output of the function eval_GASP_RFP.
}
}
\value{
Function returns a list of three objects
\item{x }{Inputs.}
\item{mu }{Mean of an emulator.}
\item{var }{Covariance matrix of an emulator.}
}
\author{
Ksenia N. Kyzyurova, kseniak.ucoz.net.
}
\examples{
## Function f2 is a simulator
f2<-function(x){cos(5*x)}

## One-dimensional inputs x2
x2 = seq(-0.95,0.95,length = 6)
data.f2 <- list(training = x2,fD = f2(x2), smooth = 2)

## Evaluation of GASP parameters
f2_MLEs = eval_GASP_RFP(data.f2,list(function(x){x^0},function(x){x^1}),1,FALSE)

## Evaluation of a T-GASP emulator
TGASP_f2 <- eval_TGASP(as.matrix(seq(-1,1,.01)),f2_MLEs)
}
\keyword{ models }
