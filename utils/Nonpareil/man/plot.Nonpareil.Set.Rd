% Hand-written Rd documentation for plot.Nonpareil.Set
\name{plot.Nonpareil.Set}
\alias{plot.Nonpareil.Set}
\title{Plot a \code{Nonpareil.Set} object.}
\usage{
\method{plot}{Nonpareil.Set}(
  x,
  col = NA,
  labels = NA,
  main = "Nonpareil Curves",
  legend.opts = list(),
  aggregate = FALSE,
  group = NULL,
  aggregate.col = NA,
  aggregate.ribbon = "sd",
  plot.individual = !aggregate,
  individual.alpha = 0.3,
  ...
)
}
\arguments{
\item{x}{\code{Nonpareil.Set} object to plot.}

\item{col}{Color of the curves (vector).
If passed, it overrides the colors set in the \code{Nonpareil.Curve}
objects. Values are recycled.}

\item{labels}{Labels of the curves (vector). If passed, it overrides the labels set
in the \code{Nonpareil.Curve} objects. Values are recycled.}

\item{main}{Title of the plot.}

\item{legend.opts}{Any additional parameters passed to \code{Nonpareil.legend}.
If FALSE, the legend is not displayed.}

\item{aggregate}{If \code{TRUE}, curves are aggregated by \code{group} using
\code{\link{Nonpareil.group}} and a consensus line (plus optional variance
ribbon) is drawn for each group instead of individual curves.  Requires
\code{group} to be supplied.  Default: \code{FALSE} (original behaviour
unchanged).}

\item{group}{Character (or coercible) vector of group labels, one per curve
in \code{x}.  Required when \code{aggregate = TRUE}.}

\item{aggregate.col}{Colors for each group (vector, recycled).  \code{NA}
entries receive random colors.  Only used when \code{aggregate = TRUE}.}

\item{aggregate.ribbon}{Type of variance ribbon drawn around the consensus
line when \code{aggregate = TRUE}.  Supported values:
\itemize{
  \item \code{"sd"} (default): ± one within-group standard deviation,
  \item \code{"ci95"}: 95\% confidence interval,
  \item \code{"ci90"}: 90\% confidence interval,
  \item \code{"ci50"}: 50\% confidence interval,
  \item \code{FALSE}: no ribbon.
}}

\item{plot.individual}{If \code{TRUE}, individual curves are drawn at
reduced alpha before the consensus overlay.  Defaults to \code{TRUE}
when \code{aggregate = FALSE} (standard behaviour) and \code{FALSE}
when \code{aggregate = TRUE}.}

\item{individual.alpha}{Alpha level (0–1) used for individual curves when
both \code{aggregate = TRUE} and \code{plot.individual = TRUE}.
Default: \code{0.3}.}

\item{...}{Any additional parameters passed to \code{plot.Nonpareil.Curve}.}
}
\value{
Returns invisibly a \code{Nonpareil.Set} object (same as \code{x} input).
}
\description{
Plot a \code{Nonpareil.Set} object.  When \code{aggregate = FALSE} (the
default), all individual curves are drawn as before.  When
\code{aggregate = TRUE}, \code{\link{Nonpareil.group}} is called internally
to compute per-group consensus curves that are then overlaid on the canvas
as solid lines with an optional variance ribbon.
}
\examples{
files <- system.file(
  "extdata",
  c("HumanGut.npo", "LakeLanier.npo", "IowaSoil.npo"),
  package = "Nonpareil"
)
col <- c("orange", "darkcyan", "firebrick4")
nps <- Nonpareil.set(
  files, col = col,
  plot.opts = list(plot.observed = FALSE, model.lwd = 2)
)

# Aggregate mode: two groups, SD ribbon
plot(nps, aggregate = TRUE, group = c("A", "A", "B"),
     aggregate.col = c("steelblue", "tomato"),
     aggregate.ribbon = "sd")

# Aggregate mode: show individual curves faintly underneath
plot(nps, aggregate = TRUE, group = c("A", "A", "B"),
     plot.individual = TRUE, individual.alpha = 0.2)
}
\seealso{
\code{\link{Nonpareil.group}}, \code{\link{plot.Nonpareil.Curve}}
}

