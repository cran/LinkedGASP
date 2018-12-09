\name{eval_GASP_RFP}
\alias{eval_GASP_RFP}
\title{
Evaluation of parameters of a Gaussian stochastic process emulator of a computer model.
}
\description{
This function evaluates parameters of a Gaussian stochastic process emulator of a computer model based on a few observations which are available from the simulator of a computer model.
}
\usage{
eval_GASP_RFP(data, basis, corr.cols, nugget)
}
\arguments{
  \item{data}{
list which consists of three objects: training input values (which may be multivariate, along several dimensions), corresponding output values of a simulator (scalar) and a vector of smoothness parameter(s) along each input direction.
}
  \item{basis}{
A set of functions in the mean of a Gaussian process. Typically assumed to be linear in one or several dimensions.
}
  \item{corr.cols}{
specifies which input directions must be included in the specification of a correlation function.
}
  \item{nugget}{
Parameter which accounts for possible small stochastisity in the output of a computer model. Default is FALSE.
}
}
\details{
See examples which illustrate inputs specification to the function.
}
\value{
%%  ~Describe the value returned
Function returns a list of objects, including estimates of parameters, which is subsequently may be used for construction of a GASP approximation with the estimated parameters and the data involved.
\item{delta }{Estimates of range parameters in the correlation function.}
\item{eta }{Estimates of a nugget.}
\item{sigma.sq }{Estimates of variance.}
\item{data }{Input parameter returned for convenience.}
\item{nugget }{Input parameter returned for convenience.}
\item{basis }{Input parameter returned for convenience.}
\item{corr.cols }{Input parameter returned for convenience.}

}
\references{Ksenia N. Kyzyurova, James O. Berger, and Robert L. Wolpert. Coupling computer models through linking their statistical emulators. SIAM/ASA Journal on Uncertainty Quantification, 6(3): 1151-1171, 2018

Gu, M., Wang, X., Berger, J. O. et al. (2018) Robust Gaussian stochastic process emulation. The Annals of Statistics, 46, 3038-3066.}
\author{
Ksenia N. Kyzyurova, kseniak.ucoz.net.
}
\examples{
## Function f1 is a simulator
f1<-function(x){sin(pi*x)}

## One-dimensional inputs are x1
x1 <- seq(-1,1,.37)

## data.f1 contains the list of data inputs (training) and outputs (fD) together with the assumed
## fixed smoothness of a computer model output. This corresponds to the smoothness in a product 
## power exponential correlation function used for construction of the emulator.
data.f1 <- list(training = x1,fD = f1(x1), smooth = 1.99)

## Evaluation of GASP parameters
f1_MLEs = eval_GASP_RFP(data.f1,list(function(x){x^0},function(x){x^1}),1,FALSE)
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
