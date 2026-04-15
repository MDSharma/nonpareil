% Hand-written Rd documentation for Nonpareil.group
\name{Nonpareil.group}
\alias{Nonpareil.group}
\title{Compute group-level consensus Nonpareil curves from a
  \code{Nonpareil.Set}.}
\usage{
Nonpareil.group(
  nps,
  group,
  col = NA,
  n.points = 500L,
  weights = c("inverse_sd", "equal"),
  extrapolation.penalty = 0.5
)
}
\arguments{
\item{nps}{\code{Nonpareil.Set} object containing the individual curves to
  aggregate.}

\item{group}{Vector of group labels (character, numeric, or factor), one
  element per curve in \code{nps}.  Curves sharing the same label are
  aggregated together.}

\item{col}{Colors for each group (vector, recycled).  \code{NA} entries are
  replaced by a random color.}

\item{n.points}{Number of points on the shared log-spaced sequencing-depth
  grid per group.  Default: \code{500L}.}

\item{weights}{Weighting scheme used when averaging across replicates at
  each depth:
  \itemize{
    \item \code{"inverse_sd"} (default): weight = 1 / within-sample SD from
      the .npo file.  Observed points carry their actual SD; projected
      (model-extrapolated) points are penalised via
      \code{extrapolation.penalty}.
    \item \code{"equal"}: all contributing curves receive weight = 1.
  }}

\item{extrapolation.penalty}{Scalar in (0, 1] applied to the weight of
  model-projected coverage values relative to observed/interpolated ones.
  A value of 0.5 (the default) halves the weight assigned to projected
  points.  Only relevant when \code{weights = "inverse_sd"}.}
}
\value{
Returns invisibly a \code{data.frame} of subclass
\code{Nonpareil.Group} with one row per (group, depth) combination and
the following columns:
\itemize{
  \item \code{group}: group label.
  \item \code{depth_bp}: sequencing depth in base pairs.
  \item \code{mean_cov}: weighted mean coverage across replicates at this
    depth.
  \item \code{sd_cov}: reliability-weighted within-group standard deviation.
    Note: uncertainty is larger where depth exceeds the observed range of
    some replicates and projection is required.
  \item \code{se_cov}: standard error of the mean.
  \item \code{ci_low}, \code{ci_high}: lower and upper bounds of the
    approximate 95\% confidence interval (mean ± 1.96 × SE).
  \item \code{n_samples}: number of curves contributing at this depth.
  \item \code{n_observed}: curves contributing interpolated (observed) data.
  \item \code{n_projected}: curves contributing model-projected data
    (coverage values beyond the observed sequencing depth of that
    replicate).
}
Object attributes \code{group.levels}, \code{group.col}, \code{nps}, and
\code{group} carry metadata used internally by \code{plot.Nonpareil.Set}
and available for custom downstream use.
}
\description{
Aggregates individual Nonpareil curves in \emph{prediction space} (not
parameter space) across replicates within each group.

For each group a log-spaced sequencing-depth grid is built spanning the
union of observed depths of its member curves.  At every grid point each
curve contributes either an interpolated value (depth within its observed
range) or a model-projected value (depth beyond observed range, requires
a fitted model).  A weighted mean and variance are then computed across
contributing curves.

The returned tidy \code{data.frame} is suitable for direct use with
\pkg{ggplot2}: map \code{depth_bp} → x, \code{mean_cov} → y, and
\code{ci_low}/\code{ci_high} to a ribbon aesthetic.

\strong{No re-running of Nonpareil is required.}  All information needed
for cross-replicate aggregation (per-depth observed coverage, within-sample
standard deviation, total library size, and model projection) is already
present in the .npo files produced by the Nonpareil binary.
}
\section{Mathematical note}{
Averaging is performed in prediction space, \emph{not} by averaging
Gamma-CDF model coefficients.  For a set of curves with parameters
\eqn{(a_i, b_i)}, the consensus coverage at depth \eqn{d} is the weighted
mean of \eqn{F(d; a_i, b_i)} values, \emph{not}
\eqn{F(d; \bar{a}, \bar{b})}.  These differ because the Gamma CDF is
nonlinear in its parameters.
}
\examples{
files <- system.file(
  "extdata",
  c("HumanGut.npo", "LakeLanier.npo", "IowaSoil.npo"),
  package = "Nonpareil"
)
nps <- Nonpareil.set(files, plot = FALSE)

# Assign two groups and compute consensus curves
grp <- Nonpareil.group(nps, group = c("A", "A", "B"))
head(grp)

# Columns available for ggplot2 ribbon plots
names(grp)

# Use the plot.Nonpareil.Set convenience wrapper
plot(nps, aggregate = TRUE, group = c("A", "A", "B"),
     aggregate.ribbon = "ci95")
}
\seealso{
\code{\link{Nonpareil.set}}, \code{\link{plot.Nonpareil.Set}},
\code{\link{predict.Nonpareil.Curve}}
}
