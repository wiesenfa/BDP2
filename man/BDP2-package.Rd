\name{BDP2-package}
\alias{BDP2-package}
\docType{package}
\title{
\packageTitle{BDP2}
}
\description{
\packageDescription{BDP2}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BDP2}
\packageIndices{BDP2}
%Workflow to choose design parameters in a adaptive Bayesian phase II trial design in a single-arm trial with binary
%endpoint (response, success) and interim efficacy and futility analyses as well as routines to determine and visualize operating characteristics.
Tools and workflow to choose design parameters in Bayesian adaptive single-arm phase II trial designs with binary
endpoint (response, success) with possible stopping for efficacy and futility at interim analyses.
Also contains routines to determine and visualize operating characteristics.


Main functions: \code{\link{BDP2workflow}}, \code{\link{BDP2}} and \code{\link{plotBDP2}}. \code{\link{BDP2workflow}} provides an interactive shiny app which also generates Word/pdf/html reports.
}
\author{
\packageAuthor{BDP2}

Maintainer: \packageMaintainer{BDP2}
}
\references{
Kopp-Schneider, A., Wiesenfarth, M., Witt, R., Edelmann, D., Witt, O. and Abel, U. (2018).\cr Monitoring futility and efficacy in phase II trials with Bayesian
posterior distributions - a calibration approach. 
\emph{Biometrical Journal}, to appear.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
\dontrun{
# Starts Shiny app
  BDP2workflow()
}
}
