\name{BDP2workflow}
\alias{BDP2workflow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shiny app for workflow
}
\description{
Starts a shiny app in the web browser. It provides a workflow to choose design parameters single-arm trial with a binary
endpoint (response, success) and interim efficacy and futility analyses as well as routines to determine and visualize operating characteristics. Also  Word/pdf/html reports can be generated.
}
\usage{
BDP2workflow(display.mode = "normal")
}
\arguments{
  \item{display.mode}{display.mode passed to \code{shiny::runApp()}}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Kopp-Schneider, A., Wiesenfarth, M., Witt, R., Edelmann, D., Witt, O. and Abel, U. (2018).\cr Monitoring futility and efficacy in phase II trials with Bayesian
posterior distributions - a calibration approach. 
\emph{Biometrical Journal}, to appear.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# Starts Shiny app
  BDP2workflow()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
