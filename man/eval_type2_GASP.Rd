\name{eval_type2_GASP}
\alias{eval_type2_GASP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The second type of an emulator of a computer model
}
\description{
This function evaluates the second GASP of a computer model within partial objective Bayesian (POB) implementation of the GASP.
}
\usage{
eval_type2_GASP(input, GASPparams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
input values (the same dimension as training input data in the next argument GASPparams)
}
  \item{GASPparams}{
The output of the function eval_GASP_RFP.
}
}
\details{
See examples which illustrate inputs specification to the function.
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

## Evaluation of a second type GASP emulator
GASP_type2_f2 <- eval_type2_GASP(as.matrix(seq(-1,1,.01)),f2_MLEs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }