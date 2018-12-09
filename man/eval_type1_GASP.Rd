\name{eval_type1_GASP}
\alias{eval_type1_GASP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The first type of an emulator of a computer model
}
\description{
This function evaluates the first GASP of a computer model using maximum a posteriori estimates (MAP) of parameters of the GASP.
}
\usage{
eval_type1_GASP(input, GASPparams)
}
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
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
