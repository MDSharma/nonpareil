# Package Docs -----------------------------------------------------------------

#' Nonpareil: Metagenome Coverage Estimation and Projections for 'Nonpareil'.
#'
#' Plot, process, and analyze NPO files produced by 'Nonpareil'
#' \url{http://enve-omics.ce.gatech.edu/nonpareil/}.
#'
#' @section Citation:
#' If you use Nonpareil, please cite:
#' Rodriguez-R et al. 2018. Nonpareil 3: Fast estimation of metagenomic coverage
#' and sequence diversity.
#' mSystems 3(3): e00039-18. DOI: 10.1128/mSystems.00039-18.
#'
#' Rodriguez-R & Konstantinidis. 2014. Nonpareil: a redundancy-based approach to
#' assess the level of coverage in metagenomic datasets.
#' Bioinformatics 30 (5): 629-635. DOI: 10.1093/bioinformatics/btt584.
#'
#' For an extended discussion on coverage in metagenomic data, see also:
#'
#' Rodriguez-R & Konstantinidis. 2014. Estimating coverage in metagenomic data
#' sets and why it matters.
#' The ISME Journal 8: 2349–2351. DOI: 10.1038/ismej.2014.76.
#'
#' @docType package
#' @name Nonpareil
NULL

# Classes ----------------------------------------------------------------------

#' A single Nonpareil curve. This object can be produced by
#' \code{Nonpareil.curve} and supports S4 methods \code{plot}, \code{summary},
#' \code{print}, and \code{predict}. For additional details, see help for
#' \code{summary.Nonpareil.Curve}.
setClass("Nonpareil.Curve",
  representation(
  
  # Dataset info
  #' @slot file Input .npo file.
  file = 'character',
  #' @slot label Name of the dataset.
  label = 'character',
  #' @slot col Color of the dataset.
  col = 'character',
  
  # Input .npo metadata
  #' @slot L Read length.
  L = 'numeric',
  #' @slot AL Adjusted read length (same as L for alignment).
  AL = 'numeric',
  #' @slot R Number of reads.
  R = 'numeric',
  #' @slot LR Effective sequencing effort used.
  LR = 'numeric',
  #' @slot overlap Minimum read overlap.
  overlap = 'numeric',
  #' @slot ksize K-mer size (for kmer kernel only).
  ksize = 'numeric',
  #' @slot log.sample Multiplier of the log-sampling (or zero if linear).
  log.sample = 'numeric',
  #' @slot kernel Read-comparison kernel.
  kernel = 'character',
  #' @slot version Nonpareil version used.
  version = 'character',

  # Input .npo data
  #' @slot x.obs Rarefied sequencing effort.
  x.obs = 'numeric',
  #' @slot x.adj Adjusted rarefied sequencing effort.
  x.adj = 'numeric',
  #' @slot y.red Rarefied redundancy (observed).
  y.red = 'numeric',
  #' @slot y.cov Rarefied coverage (corrected).
  y.cov = 'numeric',
  #' @slot y.sd Standard deviation of rarefied coverage.
  y.sd = 'numeric',
  #' @slot y.p25 Percentile 25 (1st quartile) of rarefied coverage.
  y.p25 = 'numeric',
  #' @slot y.p50 Percentile 50 (median) of rarefied coverage.
  y.p50 = 'numeric',
  #' @slot y.p75 Percentile 75 (3rd quartile) of rarefied coverage.
  y.p75 = 'numeric',

  # Estimated coverage
  #' @slot kappa Dataset redundancy.
  kappa = 'numeric',
  #' @slot C Dataset coverage.
  C = 'numeric',
  #' @slot consistent Is the data sufficient for accurate estimation?
  consistent = 'logical',

  # Projected coverage
  #' @slot star Coverage considered 'nearly complete'.
  star = 'numeric',
  #' @slot has.model Was the model successfully estimated?
  has.model = 'logical',
  #' @slot warning Warnings generated on consistency or model fitting.
  warning = 'character',
  #' @slot LRstar Projected seq. effort for nearly complete coverage.
  LRstar = 'numeric',
  #' @slot modelR Pearson's R for the estimated model.
  modelR = 'numeric',
  #' @slot diversity Dataset Nd index of sequence diversity.
  diversity = 'numeric',
  #' @slot model Fitted sigmoidal model.
  model = 'list',
  # Call
  #' @slot call Call producing this object.
  call = 'call'

  ), package = "Nonpareil"
);

#' Collection of \code{Nonpareil.Curve} objects. This object can be produced by
#' \code{Nonpareil.curve.batch} and supports S4 methods \code{plot},
#' \code{summary}, and \code{print}.
setClass(
  "Nonpareil.Set",
  representation(
    #' @slot np.curves List of \code{Nonpareil.Curve} objects.
    np.curves = "list",
    #' @slot call Call producing this object.
    call = "call"

  ), package = "Nonpareil"
);

# S3 Methods -------------------------------------------------------------------

#' Get attribute.
setMethod("$", "Nonpareil.Curve",
  #' @param x \code{Nonpareil.Curve} object.
  #' @param name Attribute.
  function(x, name) attr(x, name))

#' Set attribute.
setMethod("$<-", "Nonpareil.Curve",
  #' @param x \code{Nonpareil.Curve} object.
  #' @param name Attribute.
  #' @param value New value.
  function(x, name, value) { attr(x, name) <- value ; x })

#' Get attribute.
setMethod("$", "Nonpareil.Set",
  #' @param x \code{Nonpareil.Set} object.
  #' @param name Attribute.
  function(x, name) attr(x, name))

#' Set attribute.
setMethod("$<-", "Nonpareil.Set",
  #' @param x \code{Nonpareil.Set} object.
  #' @param name Attribute.
  #' @param value New value.
  function(x, name, value) { attr(x, name) <- value ; x })

#' Plot a \code{Nonpareil.Set} object.
plot.Nonpareil.Set <- function(
      #' @param x
      #' \code{Nonpareil.Set} object to plot.
      x,
      #' @param col
      #' Color of the curves (vector).
      #' If passed, it overrides the colors set in the \code{Nonpareil.Curve}
      #' objects. Values are recycled.
      col = NA,
      #' @param labels
      #' Labels of the curves (vector). If passed, it overrides the labels set
      #' in the \code{Nonpareil.Curve} objects. Values are recycled.
      labels = NA,
      #' @param main
      #' Title of the plot.
      main = "Nonpareil Curves",
      #' @param legend.opts
      #' Any additional parameters passed to \code{Nonpareil.legend}.
      #' If FALSE, the legend is not displayed.
      legend.opts = list(),
      #' @param aggregate
      #' If \code{TRUE}, curves are aggregated by \code{group} using
      #' \code{Nonpareil.group} and a consensus line (plus optional variance
      #' ribbon) is drawn for each group instead of individual curves.
      #' Requires \code{group} to be supplied.  Default: \code{FALSE}.
      aggregate = FALSE,
      #' @param group
      #' Character (or coercible) vector of group labels, one per curve in
      #' \code{x}.  Required when \code{aggregate = TRUE}.
      group = NULL,
      #' @param aggregate.col
      #' Colors for each group (vector).  Recycled; \code{NA} assigns random
      #' colors.  Only used when \code{aggregate = TRUE}.
      aggregate.col = NA,
      #' @param aggregate.ribbon
      #' Type of variance ribbon drawn around the consensus line when
      #' \code{aggregate = TRUE}.  Supported values:
      #' \itemize{
      #'   \item \code{"sd"} (default): ± one within-group SD,
      #'   \item \code{"ci95"}: 95\% confidence interval,
      #'   \item \code{"ci90"}: 90\% confidence interval,
      #'   \item \code{"ci50"}: 50\% confidence interval,
      #'   \item \code{FALSE}: no ribbon.
      #' }
      aggregate.ribbon = "sd",
      #' @param plot.individual
      #' If \code{TRUE}, individual curves are drawn at reduced alpha before
      #' the consensus overlay.  Defaults to \code{TRUE} when
      #' \code{aggregate = FALSE} (standard behaviour) and \code{FALSE} when
      #' \code{aggregate = TRUE}.
      plot.individual = !aggregate,
      #' @param individual.alpha
      #' Alpha level (0–1) used for individual curves when both
      #' \code{aggregate = TRUE} and \code{plot.individual = TRUE}.
      #' Default: \code{0.3}.
      individual.alpha = 0.3,
      #' @param ...
      #' Any additional parameters passed to \code{plot.Nonpareil.Curve}.
      ...
      ) {
  if (!inherits(x, "Nonpareil.Set"))
    stop("'x' must inherit from class `Nonpareil.Set`")

  # Apply per-curve col / label overrides
  col    <- rep(col,    length.out = length(x$np.curves))
  labels <- rep(labels, length.out = length(x$np.curves))
  for (i in seq_along(x$np.curves)) {
    if (!is.na(col[i]))    x$np.curves[[i]]$col   <- col[i]
    if (!is.na(labels[i])) x$np.curves[[i]]$label <- labels[i]
  }

  if (!aggregate) {
    # --- Original behaviour (unchanged) ---
    new <- TRUE
    for (i in seq_along(x$np.curves)) {
      plot(x$np.curves[[i]], new = new, main = ifelse(new, main, ""), ...)
      new <- FALSE
    }
    if (inherits(legend.opts, "list")) {
      legend.opts[["np"]] <- x
      if (is.null(legend.opts[["x"]])) legend.opts[["x"]] <- "bottomright"
      do.call(Nonpareil.legend, legend.opts)
    }
    return(invisible(x))
  }

  # --- Aggregate mode ---
  if (is.null(group))
    stop("'group' must be provided when 'aggregate = TRUE'")
  if (length(group) != length(x$np.curves))
    stop("'group' must have the same length as the number of curves in 'x'")

  grp_data     <- Nonpareil.group(x, group, col = aggregate.col)
  group.levels <- attr(grp_data, "group.levels")
  group.col    <- attr(grp_data, "group.col")

  dots <- list(...)

  if (plot.individual) {
    # Draw individual curves at reduced alpha for context
    new <- TRUE
    for (i in seq_along(x$np.curves)) {
      i_args <- c(
        list(x     = x$np.curves[[i]],
             new   = new,
             main  = ifelse(new, main, ""),
             curve.alpha = individual.alpha,
             model.alpha = individual.alpha),
        dots[!names(dots) %in% c("curve.alpha", "model.alpha")]
      )
      do.call(plot, i_args)
      new <- FALSE
    }
  } else {
    # Initialise canvas only — suppress data drawing
    canvas_args <- c(
      list(x               = x$np.curves[[1]],
           new             = TRUE,
           main            = main,
           plot.observed   = FALSE,
           plot.model      = FALSE,
           plot.diversity  = FALSE),
      dots[!names(dots) %in%
             c("plot.observed", "plot.model", "plot.diversity", "main")]
    )
    do.call(plot, canvas_args)
  }

  # Determine ylim for ribbon clamping
  ylim <- if (!is.null(dots$ylim)) dots$ylim else c(1e-6, 1)

  # Overlay ribbons then consensus lines
  for (grp in group.levels) {
    gd <- grp_data[grp_data$group == grp & !is.na(grp_data$mean_cov), ]
    if (nrow(gd) == 0L) next
    gcol <- group.col[[grp]]

    # Ribbon
    if (!isFALSE(aggregate.ribbon)) {
      lo <- hi <- NULL
      if (aggregate.ribbon == "sd") {
        lo <- gd$mean_cov - gd$sd_cov
        hi <- gd$mean_cov + gd$sd_cov
      } else if (aggregate.ribbon == "ci95") {
        lo <- gd$ci_low
        hi <- gd$ci_high
      } else if (aggregate.ribbon == "ci90") {
        lo <- gd$mean_cov - 1.645 * gd$se_cov
        hi <- gd$mean_cov + 1.645 * gd$se_cov
      } else if (aggregate.ribbon == "ci50") {
        lo <- gd$mean_cov - 0.674 * gd$se_cov
        hi <- gd$mean_cov + 0.674 * gd$se_cov
      }
      if (!is.null(lo) && !is.null(hi)) {
        lo <- pmax(lo, ylim[1] * 0.1, na.rm = TRUE)
        hi <- pmin(hi, ylim[2],        na.rm = TRUE)
        ribbon.col <- apply(
          col2rgb(gcol), 2L,
          function(ch) do.call(rgb, as.list(c(ch[1:3] / 256, 0.2)))
        )
        polygon(
          c(gd$depth_bp, rev(gd$depth_bp)),
          c(hi, rev(lo)),
          col = ribbon.col, border = NA
        )
      }
    }

    # Consensus line
    lines(gd$depth_bp, gd$mean_cov, col = gcol, lwd = 2)
  }

  # Legend for groups
  if (inherits(legend.opts, "list")) {
    leg.x <- if (!is.null(legend.opts[["x"]])) legend.opts[["x"]] else "bottomright"
    leg.y <- if (!is.null(legend.opts[["y"]])) legend.opts[["y"]] else 0.3
    extra.opts <- legend.opts[!names(legend.opts) %in% c("np", "x", "y")]
    do.call(legend, c(
      list(x = leg.x, y = leg.y,
           legend = group.levels,
           fill   = unname(group.col)),
      extra.opts
    ))
  }

  #' @return
  #' Returns invisibly a \code{Nonpareil.Set} object (same as \code{x} input).
  invisible(x)
}

#' Plot a \code{Nonpareil.Curve} object.
plot.Nonpareil.Curve <- function(
      #' @param x
      #' \code{Nonpareil.Curve} object to plot.
      x,
      #' @param col
      #' Color of the curve. If passed, it overrides the colors set in the
      #' \code{Nonpareil.Curve} object.
      col = NA,
      #' @param add
      #' If TRUE, it attempts to use a previous (active) canvas to plot the
      #' curve.
      add = FALSE,
      #' @param new
      #' Inverse of `add`.
      new = !add,
      #' @param plot.observed
      #' Indicates if the observed (rarefied) coverage is to be plotted.
      plot.observed = TRUE,
      #' @param plot.model
      #' Indicates if the fitted model is to be plotted.
      plot.model = TRUE,
      #' @param plot.dispersion
      #' Indicates if (and how) dispersion of the replicates should be
      #' plotted. Supported values are:
      #' \itemize{
      #'   \item \code{FALSE}: no dispersion is plotted (default),
      #'   \item \code{'sd'}: one standard deviation around the mean,
      #'   \item \code{'ci95'}: 95% confidence interval,
      #'   \item \code{'ci90'}: 90% confidence interval,
      #'   \item \code{'ci50'}: 50% confidence interval,
      #'   \item \code{'iq'}: Inter-quartile range.
      #' }
      plot.dispersion = FALSE,
      #' @param plot.diversity
      #' If TRUE, the diversity estimate is plotted as a small arrow below the
      #' Nonpareil curve.
      plot.diversity = TRUE,
      #' @param xlim
      #' Limits of the sequencing effort (X-axis).
      xlim = c(1e3, 10e12),
      #' @param ylim
      #' Limits of the coverage (Y-axis).
      ylim = c(1e-6, 1),
      #' @param main
      #' Title of the plot.
      main = paste("Nonpareil Curve for", x$label),
      #' @param xlab
      #' Label of the X-axis.
      xlab = "Sequencing effort (bp)",
      #' @param ylab
      #' Label of the Y-axis.
      ylab = "Estimated Average Coverage",
      #' @param curve.lwd
      #' Line width of the rarefied coverage.
      curve.lwd = 2,
      #' @param curve.alpha
      #' Alpha value (from 0 to 1) of the rarefied coverage.
      curve.alpha=0.4,
      #' @param model.lwd
      #' Line width of the model.
      model.lwd = 1,
      #' @param model.alpha
      #' Alpha value (from 0 to 1) of the model.
      model.alpha = 1,
      #' @param log
      #' Axis to plot in logarithmic scale. Supported values are:
      #' \itemize{
      #'   \item \code{'x'}: sequencing effort (default),
      #'   \item \code{'y'}: coverage,
      #'   \item \code{'xy'}: both logarithmic, or
      #'   \item \code{''}: both linear.
      #' }
      log = "x",
      #' @param arrow.length
      #' If \code{plot.diversity = TRUE}, it determines the length of the
      #' arrow to display the divsersity (as a fraction of the ylim range).
      arrow.length = 0.05,
      #' @param arrow.head
      #' If \code{plot.diversity = TRUE}, it determines the length of the
      #' arrow head to display the diversity index (in inches).
      arrow.head = arrow.length,
      #' @param ...
      #' Additional graphical parameters.
      ...
      ){
  if(!inherits(x, "Nonpareil.Curve"))
    stop("'x' must inherit from class `Nonpareil.Curve`")

  # Create empty canvas
  if(new){
    plot(
      1, type = "n", xlim = xlim, ylim = ylim, bty = "l",
      xlab = xlab, ylab = ylab, main = main,
      xaxs = "i", yaxs = "i", log = log, ...
    )
    abline(h = c(1, x$star/100), lty = 2, col = "red")
    abline(v = 10^seq(0, 15, by = 3), lty = 2, col = "gray80")
  }

  # Dispersion
  if(plot.dispersion!=FALSE){
    if(plot.dispersion == "sd"){
      err.y <- c(x$y.cov+x$y.sd, rev(x$y.cov-x$y.sd))
    }else if(plot.dispersion == "ci95"){
      err.y <- c(x$y.cov+x$y.sd*1.9, rev(x$y.cov-x$y.sd*1.9))
    }else if(plot.dispersion == "ci90"){
      err.y <- c(x$y.cov+x$y.sd*1.64, rev(x$y.cov-x$y.sd*1.64))
    }else if(plot.dispersion == "ci50"){
      err.y <- c(x$y.cov+x$y.sd*.67, rev(x$y.cov-x$y.sd*.67))
    }else if(plot.dispersion == "iq"){
      err.y <- c(x$y.p25, rev(x$y.p75))
    }
    polygon(c(x$x.adj, rev(x$x.adj)), border=NA,
      ifelse(err.y<=ylim[1]*0.1, ylim[1]*0.1, err.y), col=Nonpareil.col(x, .2))
  }

  # Rarefied coverage
  if(plot.observed){
    lines(x$x.adj, x$y.cov, col=Nonpareil.col(x, curve.alpha), lwd=curve.lwd);
  }

  # Model
  if(x$has.model & plot.model){
    model.lty <- ifelse(plot.observed, 2, 1)
    model.x   <- exp(seq(log(xlim[1]), log(xlim[2]), length.out=1e3));
    model.y   <- predict(x, lr=model.x);
    lines(model.x, model.y, col=Nonpareil.col(x, model.alpha), lty=model.lty,
          lwd=model.lwd)
    if(!plot.observed){
      points(x$LR, predict(x), col=Nonpareil.col(x, 1.0), pch=21, bg="white")
    }
  }

  if(x$has.model & plot.diversity & x$diversity>0){
    arrows(x0=exp(x$diversity), length=arrow.head,
          y1=ifelse(log=='y' | log=='xy' | log=='yx',
            ylim[1]*(ylim[2]/ylim[1])**arrow.length,
            ylim[1] + diff(ylim)*arrow.length),
          y0=ylim[1], col=Nonpareil.col(x, model.alpha));
  }

  #' @return
  #' Retuns invisibly a \code{Nonpareil.Curve} object (same as \code{x} input).
  #' For additional details see help for \code{summary.Nonpareil.Curve}.
  invisible(x)
}

#' Returns a summary of the \code{Nonpareil.Set} results.
summary.Nonpareil.Set <- function(
      #' @param object
      #' \code{Nonpareil.Set} object.
      object,
      #' @param ...
      #' Additional parameters ignored.
      ...
      ){
  if(!inherits(object, "Nonpareil.Set"))
    stop("'object' must inherit from class `Nonpareil.Set`")
  y <- rbind(sapply(object$np.curves, "summary"))
  colnames(y) <- sapply(object$np.curves, function(n) n$label)

  #' @return
  #' Returns a matrix with different values for each dataset. For additional
  #' details on the values returned, see help for
  #' \code{summary.Nonpareil.Curve}.
  t(y)
}

#' Returns a summary of the \code{Nonpareil.Curve} results.
summary.Nonpareil.Curve <- function(
      #' @param object
      #' \code{Nonpareil.Curve} object.
      object,
      #' @param ...
      #' Additional parameters ignored.
      ...
      ){
  if(!inherits(object, "Nonpareil.Curve"))
    stop("'object' must inherit from class `Nonpareil.Curve`")
  n <- c("kappa","C","LR","modelR","LRstar","diversity")
  y <- sapply(n, function(v) attr(object,v))
  names(y) <- n

  #' @return
  #' Returns a matrix with the following values for the dataset:
  #' \itemize{
  #'   \item kappa: "Redundancy" value of the entire dataset.
  #'   \item C: Average coverage of the entire dataset.
  #'   \item LRstar: Estimated sequencing effort required to reach the objective
  #'     average coverage (star, 95% by default).
  #'   \item LR: Actual sequencing effort of the dataset.
  #'   \item modelR: Pearson's R coefficient betweeen the rarefied data and the
  #'     projected model.
  #'   \item diversity: Nonpareil sequence-diversity index (Nd). This value's
  #'     units are the natural logarithm of the units of sequencing effort
  #'     (log-bp), and indicates the inflection point of the fitted model for
  #'     the Nonpareil curve. If the fit doesn't converge, or the model is not
  #'     estimated, the value is zero (0).
  #' }
  y
}

#' Prints and returns invisibly a summary of the \code{Nonpareil.Set} results.
print.Nonpareil.Set <- function(
      #' @param x
      #' \code{Nonpareil.Set} object.
      x,
      #' @param ...
      #' Additional parameters ignored.
      ...
      ){
  if(!inherits(x, "Nonpareil.Set"))
    stop("'x' must inherit from class `Nonpareil.Set`")
  y <- summary(x)
  cat("===[ Nonpareil.Set ]===================================\n")
  cat("Collection of", length(x$np.curves), "Nonpareil curves.\n")
  print(y)
  cat("-------------------------------------------------------\n")
  cat("call:",as.character(x$call),"\n")
  cat("-------------------------------------------------------\n")

  #' @return
  #' Returns the summary invisibly. See help for
  #' \code{summary.Nonpareil.Curve} and \code{summary.Nonpareil.Set} for
  #' additional information.
  invisible(y)
}

#' Prints and returns invisibly a summary of the \code{Nonpareil.Curve} results.
print.Nonpareil.Curve <- function(
      #' @param x
      #' \code{Nonpareil.Set} object.
      x,
      #' @param ...
      #' Additional parameters ignored.
      ...
      ){
  if(!inherits(x, "Nonpareil.Curve"))
    stop("'x' must inherit from class `Nonpareil.Curve`")
  y <- summary(x)
  yp <- cbind(y)
  colnames(yp) <- x$label
  cat("===[ Nonpareil.Curve ]=================================\n")
  print(yp)
  cat("-------------------------------------------------------\n")
  cat("call:",as.character(x$call),"\n")
  cat("-------------------------------------------------------\n")

  #' @return
  #' Returns the summary invisibly. See help for
  #' \code{summary.Nonpareil.Curve} for additional information.
  invisible(y)
}

#' Predict the coverage for a given sequencing effort.
predict.Nonpareil.Curve <- function(
      #' @param object
      #' \code{Nonpareil.Curve} object.
      object,
      #' @param lr
      #' Sequencing effort for the prediction (in bp).
      lr = object$LR,
      #' @param ...
      #' Additional parameters ignored.
      ...
      ){
  if(!inherits(object, "Nonpareil.Curve"))
    stop("'object' must inherit from class `Nonpareil.Curve`")
  if(!object$has.model)
    stop("'object' must be a Nonpareil Curve with a fitted model")
  
  #' @return
  #' Returns the expected coverage at the given sequencing effort.
  predict(object$model, list(x = lr))
}

# Ancillary functions ----------------------------------------------------------

#' Read the metadata headers.
Nonpareil.read_metadata <- function(
      #' @param x
      #' \code{Nonpareil.Curve} object.
      x
      ){
  # Load key-values and defaults
  meta_data <- gsub('^# @', "", grep("^# @", readLines(x$file), value = TRUE))
  keys <- gsub(': .*', "", meta_data)
  vals <- gsub('.*: ', "", meta_data)
  x$kernel <- "alignment"
  x$log.sample <- 0

  # Set metadata
  if("ksize" %in% keys && vals[keys=="ksize"]>0 && vals[keys=="ksize"]<1001)
    x$kernel <- "kmer";
  if("divide" %in% keys)
    x$log.sample <- as.numeric(vals[keys=="divide"]);
  if("logsampling" %in% keys)
    x$log.sample <- as.numeric(vals[keys=="logsampling"]);
  x$version   <- as.character(vals[keys=="version"])
  x$L         <- as.numeric(vals[keys=="L"])
  x$R         <- as.numeric(vals[keys=="R"])
  if(x$kernel=="kmer"){
    x$overlap <- 50
    x$ksize   <- as.numeric(vals[keys=="ksize"])
    x$AL      <- as.numeric(vals[keys=="AL"])
  }else{
    x$overlap <- as.numeric(vals[keys=="overlap"])
    x$AL      <- x$L
  }
  x$LR <- exp(log(x$R) + log(x$L));

  invisible(x)
}

#' Read the data tables and extract direct estimates.
Nonpareil.read_data <- function(
      #' @param x
      #' \code{Nonpareil.Curve} object.
      x,
      #' @param correction.factor
      #' Logical; see \code{Nonpareil.curve} for details.
      correction.factor
      ){
  # Read input
  a <- read.table(x$file, sep="\t", header=FALSE)
  a <- a[order(a[,1]),]
  x$x.obs <- a[,1]
  x$y.red <- a[,2]
  x$kappa <- tail(x$y.red, n=1)

  # Estimate coverage
  cor.f   <- 1.0;
  if(correction.factor) cor.f <- Nonpareil.coverage_factor(x)
  for(i in 2:6) a[, i] <- a[, i]^cor.f;
  x$y.cov <- a[, 2]
  x$y.sd  <- a[, 3]
  x$y.p25 <- a[, 4]
  x$y.p50 <- a[, 5]
  x$y.p75 <- a[, 6]
  x$C <- tail(x$y.cov, n=1)

  # Adjust sequencing effort
  x$x.adj <- exp(
    max(log(x$x.obs)) + (x$C^0.27)*(log(x$x.obs) - max(log(x$x.obs)))
  )
  # Obsolete corrections {
  #   x$x.adj <- exp(log(x$x.adj)*0.61 + 10)
  #   x$x.adj <- x$x.adj * x$AL / 101
  # }
  x$x.adj <- x$x.adj * x$AL * x$R / max(x$x.adj)

  # Check consistency
  x$consistent <- TRUE
  twenty.pc = which.max(x$x.adj[x$x.adj <= 0.5*tail(x$x.adj, n=1)]);
  if(length(twenty.pc) == 0) twenty.pc = length(x$x.adj)
  if(x$y.p50[twenty.pc] == 0){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
        paste("Median of the curve is zero at 50% of the reads, check",
        "parameters and re-run (e.g., decrease -L in nonpareil -T alignment)."))
  }
  if(x$kappa <= 1e-5){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
        paste("Redundancy curve too low, check parameters and re-run",
        "(e.g., decrease -L in nonpareil -T alignment)."))
  }
  if(x$y.cov[2] >= 1-1e-5){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
        paste("Curve too steep, check parameters and re-run",
        "(e.g., increase value of -d in nonpareil)."))
  }
  if(sum(x$y.cov>0 & x$y.cov<0.9) <= 10){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
        paste("Insufficient resolution below 90% coverage, check",
        "parameters and re-run (e.g., increase the value of -d in nonpareil)."))
  }

  invisible(x)
}

#' Fit the sigmoidal model to the rarefied coverage.
Nonpareil.fit_model <- function(
      #' @param np
      #' \code{Nonpareil.Curve} object.
      np,
      #' @param weights.exp
      #' Numeric; see \code{Nonpareil.curve} for details.
      weights.exp
      ){
  if(!inherits(np, "Nonpareil.Curve"))
    stop("'np' must inherit from class `Nonpareil.Curve`")

  # Prepare data
  sel <- np$y.cov>0 & np$y.cov<0.9
  data <- list(x=np$x.adj[ sel ], y=np$y.cov[ sel ])
  if(is.na(weights.exp[1])){
    if(np$log.sample==0){ weights.exp <- c(-1.1,-1.2,-0.9,-1.3,-1) }
    else{ weights.exp <- c(0,1,-1,1.3,-1.1,1.5,-1.5,3,-3) }
  }

  # Find the first weight with proper fit
  np$has.model <- FALSE
  weights.i <- 0
  while(!np$has.model & !is.na(weights.exp[weights.i+1])){
    weights.i <- weights.i+1
    suppressWarnings(
      model <- nls(y ~ Nonpareil.f(x, a, b), data=data,
            weights=(np$y.sd[sel]^weights.exp[weights.i]),
            start=list(a=1, b=0.1), lower=c(a=0, b=0), algorithm="port",
             control=nls.control(
                   minFactor=1e-25000, tol=1e-15, maxiter=1024, warnOnly=TRUE))
    )
    tryCatch({ is.conv <- summary(model)$convInfo$isConv },
          error=function(e){ is.conv <- FALSE })
    if(is.conv){
      np$model <- model
      np$has.model <- TRUE
    }
  }

  # Estimate diversity and projections
  if(np$has.model){
    pa <- coef(np$model)[1]
    pb <- coef(np$model)[2]
    if(pa > 1) np$diversity <- (pa-1)/pb
    np$LRstar <- Nonpareil.antif(np$star/100, pa, pb)
    np$modelR <- cor(data$y, predict(np, lr=data$x))
  }else{
    np$warnings <- c(np$warnings,
          "Model didn't converge. Try modifying the values of weights.exp.")
  }

  invisible(np)
}

#' Factor to transform redundancy into coverage (internal function).
Nonpareil.coverage_factor <- function(
      #' @param x
      #' \code{Nonpareil.Curve} object.
      x
      ){
  #' @return
  #' A numeric scalar.
  return(1 - exp(2.23E-2 * x$overlap - 3.5698))
}

#' Returns the color of the curve.
Nonpareil.col <- function(
      #' @param x
      #' \code{Nonpareil.Curve} or \code{Nonpareil.Set} object.
      x,
      #' @param alpha
      #' Alpha level of the color from 0 to 1.
      alpha = 1
      ){
  if(inherits(x, "Nonpareil.Curve")){
    col <- x$col
  }else if(inherits(x, "Nonpareil.Set")){
    col <- sapply(x$np.curves, function(np) np$col)
  }else{
    stop("'x' must inherit from class `Nonpareil.Curve` or `Nonpareil.Set`")
  }
  apply(col2rgb(col), 2,
    function(x) do.call(rgb, as.list(c(x[1:3]/256, alpha))) )
}

#' Generates a legend for Nonpareil plots.
Nonpareil.legend <- function(
      #' @param np
      #' A \code{Nonpareil.Set} or \code{Nonpareil.Curve} object, or a list of
      #' \code{Nonpareil.Curve} objects.
      np,
      #' @param x
      #' X coordinate, or any character string accepted by legend (e.g.,
      #' 'bottomright').
      x,
      #' @param y
      #' Y coordinate.
      y = 0.3,
      #' @param ...
      #' Any other parameters supported by legend().
      ...
      ){
  if(inherits(np, "Nonpareil.Set")) np <- np$np.curves
  if(inherits(np, "Nonpareil.Curve")) np <- list(np)
  if(!inherits(np, "list"))
    stop("'np' must inherit from `list`, `Nonpareil.Set`, or `Nonpareil.Curve`")
  if(missing(x)) x <- 'bottomright'

  labels <- sapply(np, function(x) x$label)
  cols <- sapply(np, Nonpareil.col)

  #' @return
  #' Returns invisibly a list, same as \code{legend}.
  legend(x = x, y = y, legend = labels, fill = cols, ...)
}

#' Adds a \code{Nonpareil.Curve} to a \code{Nonpareil.Set}.
Nonpareil.add.curve <- function(
      #' @param nps
      #' \code{Nonpareil.Set} object.
      nps,
      #' @param np
      #' \code{Nonpareil.Curve} object.
      np
      ){
  if(!inherits(nps, "Nonpareil.Set"))
    stop("'nps' must inherit from class `Nonpareil.Set`")
  if(!inherits(np, "Nonpareil.Curve"))
    stop("'np' must inherit from class `Nonpareil.Curve`")

  nps$np.curves[[ length(nps$np.curves)+1 ]] <- np

  #' @return
  #' Returns the \code{Nonpareil.Set} including the newly added
  #' \code{Nonpareil.Curve}.
  return(nps)
}

#' Alias of \code{Nonpareil.add.curve}.
setMethod("+", "Nonpareil.Set",
  #' @param e1 \code{Nonpareil.Set} object (\code{nps}).
  #' @param e2 \code{Nonpareil.Curve} object (\code{np}).
  function(e1,e2) Nonpareil.add.curve(e1,e2))

# Model Functions --------------------------------------------------------------

#' Function of the projected model.
Nonpareil.f <- function(
      #' @param x
      #' Values of sequencing effort (in bp).
      x,
      #' @param a
      #' Parameter alpha of the Gamma CDF.
      a,
      #' @param b
      #' Parameter beta of the Gamma CDF.
      b
      ){
  #' @return
  #' Predicted values of abundance-weighted average coverage.
  return(pgamma(log1p(x), a, b))
}

#' Complement function of \code{Nonpareil.f}.
Nonpareil.antif <- function(
    #' @param y
    #' Values of abundance-weighted average coverage.
    y,
    #' @param a
    #' Parameter alpha of the gamma CDF.
    a,
    #' @param b
    #' Parameter beta of the gamma CDF.
    b
    ){
  #' @return
  #' Estimated sequencing effort.
  return(exp(qgamma(y,a,b))-1)
}

# Main functions ---------------------------------------------------------------

#' Generates a Nonpareil curve from an .npo file
Nonpareil.curve <- function(
      #' @param file
      #' Path to the .npo file, containing the read redundancy.
      file,
      #' @param plot
      #' Determines if the plot should be produced. If FALSE, it computes the
      #' coverage and the model wihtout plotting.
      plot = TRUE,
      #' @param label
      #' Name of the dataset. If NA, it is determined by the file name.
      label = NA,
      #' @param col
      #' Color of the curve.
      #' If NA, a random color is assigned (even if \code{plot = FALSE}).
      col = NA,
      #' @param enforce.consistency
      #' If TRUE, it fails verbosely on insufficient data, otherwise it warns
      #' about the inconsistencies and attempts the estimations.
      enforce.consistency = TRUE,
      #' @param star
      #' Objective coverage in percentage; i.e., coverage value considered
      #' near-complete.
      star = 95,
      #' @param correction.factor
      #' Should the overlap-dependent (or kmer-length-dependent) correction
      #' factor be applied? If FALSE, redundancy is assumed to equal coverage.
      correction.factor = TRUE,
      #' @param weights.exp
      #' Vector of values to be tested (in order) as exponent of the weights
      #' distribution. If the model fails to converge, sometimes manual
      #' modifications in this parameter may help. By default (NA), five
      #' different values are tested in the following order: For linear
      #' sampling, -1.1, -1.2, -0.9, -1.3, -1. For logarithmic sampling (-d
      #' option in Nonpareil), 0, 1, -1, 1.3, -1.1, 1.5, -1.5.
      weights.exp = NA,
      #' @param skip.model
      #' If set, skips the model estimation altogether.
      skip.model = FALSE,
      #' @param ...
      #' Any additional parameters passed to \code{plot.Nonpareil.Curve}.
      ...
      ){
  # Check parameters and initialize object
  if(is.na(label)){
    label <- basename(file)
    if(substr(label, nchar(label)-3, nchar(label))==".npo"){
      label <- substr(label, 0, nchar(label)-4)
    }
  }
  if(is.na(col)){
    col <- rgb(sample(200,1),sample(200,1),sample(200,1),maxColorValue=255)
  }
  np <- new("Nonpareil.Curve", file=as.character(file), label=label, col=col,
        star=star, has.model=FALSE, call=match.call(), diversity=0)

  # Read metadata (.npo headers)
  np <- Nonpareil.read_metadata(np)

  # Read data (.npo table)
  np <- Nonpareil.read_data(np, correction.factor)
  if(!np$consistent & enforce.consistency){
    for(w in np$warnings) warning(w)
    return(invisible(np))
  }

  # Fit model
  if(!skip.model) np <- Nonpareil.fit_model(np, weights.exp)

  # Plot
  if(plot) plot(np, ...)

  # Warnings
  for(w in np$warnings) warning(w)

  #' @return
  #' Returns invisibly a \code{Nonpareil.Curve} object
  invisible(np)

  #' @examples
  #' # Generate a Nonpareil plot
  #' file <- system.file("extdata", "LakeLanier.npo", package = "Nonpareil")
  #' np <- Nonpareil.curve(file)
  #'
  #' # Produce the same plot but using powers of 1,000bp as X axis labels
  #' Nonpareil.curve(file, xaxt = "n", xlab = "Sequencing Effort")
  #' axis(
  #'   1L, at = 10L^seq(3L, 12L, by = 3L),
  #'   labels = paste(1L, c("Kbp", "Mbp", "Gbp", "Tbp"))
  #' )
  #'
  #' # Show the estimated values
  #' print(np)
  #'
  #' # Predict coverage for 20Gbp
  #' predict(np, 20e9)
  #'
  #' # Obtain the Nd diversity index
  #' np$diversity
}

#' Generates a collection of Nonpareil curves (a \code{Nonpareil.Set} object)
#' and (optionally) plots all of them in a single canvas.
Nonpareil.set <- function(
      #' @param files
      #' Vector with the paths to the .npo files.
      files,
      #' @param col
      #' Color of the curves (vector). If not passed, values are randomly
      #' assigned. Values are recycled.
      col = NA,
      #' @param labels
      #' Labels of the curves (vector). If not passed, values are determined by
      #' the filename. Values are recycled.
      labels = NA,
      #' @param plot
      #' If TRUE, it generates the Nonpareil curve plots.
      plot = TRUE,
      #' @param plot.opts
      #' Any parameters accepted by \code{plot.Nonpareil.Set} as a list.
      plot.opts = list(),
      #' @param ...
      #' Any additional parameters accepted by \code{Nonpareil.curve}.
      ...
      ){
  files <- as.vector(files)
  y <- new("Nonpareil.Set", call=match.call())
  col <- rep(col, length.out=length(files))
  labels <- rep(labels, length.out=length(files))

  nonpareil.opts <- list(...)
  nonpareil.opts[["plot"]] <- FALSE
  for(i in 1:length(files)){
    nonpareil.opts[["file"]] <- files[i]
    nonpareil.opts[["col"]] <- col[i]
    nonpareil.opts[["label"]] <- labels[i]
    y$np.curves[[ length(y$np.curve)+1 ]] <- do.call("Nonpareil.curve",
      nonpareil.opts)
  }

  # Plot
  if(plot){
    plot.opts[["x"]] <- y
    y <- do.call("plot", plot.opts)
  }

  #' @return
  #' Returns invisibly a \code{Nonpareil.Set} object.
  invisible(y)

  #' @examples
  #' # Generate a Nonpareil plot with multiple curves
  #' files <- system.file(
  #'   "extdata",
  #'   c("HumanGut.npo", "LakeLanier.npo", "IowaSoil.npo"),
  #'   package = "Nonpareil"
  #' )
  #' col <- c("orange","darkcyan","firebrick4")
  #' nps <- Nonpareil.set(
  #'   files, col = col,
  #'   plot.opts = list(plot.observed = FALSE, model.lwd = 2)
  #' )
  #'
  #' # Show the estimated values
  #' print(nps)
  #'
  #' # Show current coverage (as %)
  #' summary(nps)[, "C"] * 100
  #'
  #' # Extract Nd diversity index
  #' summary(nps)[, "diversity"]
  #'
  #' # Extract sequencing effort for nearly complete coverage (in Gbp)
  #' summary(nps)[, "LRstar"] / 1e9
  #'
  #' # Predict coverage for a sequencing effort of 10Gbp
  #' sapply(nps$np.curves, predict, 10e9)
}

#' Alias of \code{Nonpareil.set}.
#' @inheritParams Nonpareil.set
Nonpareil.curve.batch <- Nonpareil.set

#' Compute group-level consensus Nonpareil curves from a \code{Nonpareil.Set}.
#'
#' Aggregates individual Nonpareil curves in *prediction space* (not parameter
#' space) across all replicates within each group.  For each group, a
#' log-spaced sequencing-depth grid is built spanning the union of observed
#' depths across its member curves.  At every grid point each curve
#' contributes either an interpolated value (when the depth falls within its
#' observed range) or a model-projected value (when the depth is beyond the
#' observed range but a fitted model exists).  A weighted mean and variance are
#' computed across contributing curves; the inverse-variance scheme uses the
#' within-sample standard deviation stored in the .npo file, so that
#' well-resolved rarefaction points receive more weight than noisier or
#' extrapolated ones.
#'
#' The function returns a \code{data.frame} (subclass
#' \code{Nonpareil.Group}) that can be used directly with \code{ggplot2}:
#' map \code{depth_bp} → x, \code{mean_cov} → y, and
#' \code{ci_low}/\code{ci_high} → a ribbon.
#'
#' @section No changes to .npo files:
#' All information required for cross-replicate aggregation (per-depth
#' observed coverage, within-sample standard deviation, library size, and
#' model projection) is already present in the .npo outputs produced by the
#' Nonpareil binary.  No re-running of Nonpareil is needed.
Nonpareil.group <- function(
      #' @param nps
      #' \code{Nonpareil.Set} object containing the individual curves to
      #' aggregate.
      nps,
      #' @param group
      #' Vector of group labels (character, numeric, or factor), one element
      #' per curve in \code{nps}.  Curves sharing the same label are
      #' aggregated together.
      group,
      #' @param col
      #' Colors for each group (vector, recycled).  \code{NA} entries are
      #' replaced by a random color.
      col = NA,
      #' @param n.points
      #' Number of points on the shared log-spaced sequencing-depth grid per
      #' group.  Default: \code{500L}.
      n.points = 500L,
      #' @param weights
      #' Weighting scheme used when averaging across replicates at each depth:
      #' \itemize{
      #'   \item \code{"inverse_sd"} (default): weight = 1 / within-sample SD
      #'     from the .npo file.  Observed points carry their actual SD;
      #'     projected (model-extrapolated) points are penalised via
      #'     \code{extrapolation.penalty}.
      #'   \item \code{"equal"}: all contributing curves receive weight = 1.
      #' }
      weights = c("inverse_sd", "equal"),
      #' @param extrapolation.penalty
      #' Scalar in (0, 1] applied to the weight of model-projected coverage
      #' values relative to observed/interpolated ones.  A value of 0.5 (the
      #' default) halves the weight assigned to projected points.  Only
      #' relevant when \code{weights = "inverse_sd"}.
      extrapolation.penalty = 0.5
      ) {
  if (!inherits(nps, "Nonpareil.Set"))
    stop("'nps' must inherit from class `Nonpareil.Set`")
  n.curves <- length(nps$np.curves)
  if (n.curves == 0L)
    stop("'nps' contains no Nonpareil curves")
  group <- as.character(group)
  if (length(group) != n.curves)
    stop("'group' must have the same length as the number of curves in 'nps'")
  weights <- match.arg(weights)
  if (extrapolation.penalty <= 0 || extrapolation.penalty > 1)
    stop("'extrapolation.penalty' must be in the range (0, 1]")

  group.levels <- unique(group)

  # Assign group colors
  col <- rep(col, length.out = length(group.levels))
  for (gi in seq_along(group.levels)) {
    if (is.na(col[gi])) {
      col[gi] <- rgb(
        sample(200L, 1L), sample(200L, 1L), sample(200L, 1L),
        maxColorValue = 255
      )
    }
  }
  group.col <- setNames(col, group.levels)

  results <- vector("list", length(group.levels))
  names(results) <- group.levels

  for (grp in group.levels) {
    idx    <- which(group == grp)
    curves <- nps$np.curves[idx]
    n.g    <- length(curves)

    # Shared log-spaced depth grid across all curves in this group
    all.x <- unlist(lapply(curves, function(np) np$x.adj[np$x.adj > 0]))
    if (length(all.x) == 0L) {
      warning("Group '", grp, "': no positive x.adj values; skipping")
      next
    }
    x.min <- min(all.x)
    x.max <- max(all.x)
    grid  <- exp(seq(log(x.min), log(x.max), length.out = n.points))

    cov.mat <- matrix(NA_real_, nrow = n.points, ncol = n.g)
    sd.mat  <- matrix(NA_real_, nrow = n.points, ncol = n.g)
    obs.mat <- matrix(FALSE,    nrow = n.points, ncol = n.g)

    for (j in seq_len(n.g)) {
      np  <- curves[[j]]
      sel <- np$x.adj > 0
      xj  <- np$x.adj[sel]
      yj  <- np$y.cov[sel]
      sj  <- np$y.sd[sel]

      # Interpolate within observed range (approx returns NA outside range)
      obs_cov <- approx(xj, yj, xout = grid, rule = 1L)$y
      obs_sd  <- approx(xj, sj, xout = grid, rule = 1L)$y

      in.range <- !is.na(obs_cov)
      cov.mat[in.range, j] <- obs_cov[in.range]
      sd.mat[in.range, j]  <- pmax(obs_sd[in.range], 0)
      obs.mat[in.range, j] <- TRUE

      # Project outside observed range using fitted model
      if (np$has.model) {
        out.range <- !in.range
        if (any(out.range)) {
          pred <- tryCatch(
            predict(np, lr = grid[out.range]),
            error = function(e) rep(NA_real_, sum(out.range))
          )
          cov.mat[out.range, j] <- pmax(0, pmin(1, pred))
          # Inflate SD for projected points: sd_proj = min_sd / penalty
          pos.sd <- sj[sj > 0]
          min.sd <- if (length(pos.sd) > 0) min(pos.sd) else .Machine$double.eps
          sd.mat[out.range, j]  <- min.sd / extrapolation.penalty
          obs.mat[out.range, j] <- FALSE
        }
      }
    }

    # Compute weights
    w.mat <- matrix(NA_real_, nrow = n.points, ncol = n.g)
    if (weights == "inverse_sd") {
      eps   <- .Machine$double.eps
      w.mat <- 1.0 / pmax(sd.mat, eps)
    } else {
      w.mat[!is.na(cov.mat)] <- 1.0
    }

    # Per-depth statistics
    n.samp <- apply(!is.na(cov.mat), 1L, sum)
    n.obs  <- apply(obs.mat,         1L, sum)
    n.proj <- n.samp - n.obs

    mean.cov <- rep(NA_real_, n.points)
    sd.cov   <- rep(NA_real_, n.points)
    se.cov   <- rep(NA_real_, n.points)

    for (i in seq_len(n.points)) {
      valid <- which(!is.na(cov.mat[i, ]))
      if (length(valid) == 0L) next
      y    <- cov.mat[i, valid]
      w    <- w.mat[i, valid]
      wsum <- sum(w)
      mean.cov[i] <- sum(w * y) / wsum
      if (length(valid) > 1L) {
        # Reliability-weighted variance (unbiased)
        v2 <- sum(w^2)
        var.w <- sum(w * (y - mean.cov[i])^2) / (wsum - v2 / wsum)
        sd.cov[i] <- sqrt(max(var.w, 0))
        se.cov[i] <- sd.cov[i] / sqrt(length(valid))
      } else {
        sd.cov[i] <- 0
        se.cov[i] <- 0
      }
    }

    mean.cov <- pmax(0, pmin(1, mean.cov))
    ci.low   <- pmax(0,          mean.cov - 1.96 * se.cov)
    ci.high  <- pmin(1,          mean.cov + 1.96 * se.cov)

    results[[grp]] <- data.frame(
      group       = grp,
      depth_bp    = grid,
      mean_cov    = mean.cov,
      sd_cov      = sd.cov,
      se_cov      = se.cov,
      ci_low      = ci.low,
      ci_high     = ci.high,
      n_samples   = n.samp,
      n_observed  = n.obs,
      n_projected = n.proj,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  if (is.null(out) || nrow(out) == 0L) {
    warning("No data produced for any group")
    out <- data.frame(
      group = character(), depth_bp = numeric(), mean_cov = numeric(),
      sd_cov = numeric(), se_cov = numeric(),
      ci_low = numeric(), ci_high = numeric(),
      n_samples = integer(), n_observed = integer(), n_projected = integer(),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL

  attr(out, "group.levels") <- group.levels
  attr(out, "group.col")    <- group.col
  attr(out, "nps")          <- nps
  attr(out, "group")        <- group
  class(out) <- c("Nonpareil.Group", "data.frame")

  #' @return
  #' Returns invisibly a \code{data.frame} of subclass
  #' \code{Nonpareil.Group} with one row per (group, depth) combination and
  #' the following columns:
  #' \itemize{
  #'   \item \code{group}: group label.
  #'   \item \code{depth_bp}: sequencing depth in base pairs.
  #'   \item \code{mean_cov}: weighted mean coverage across replicates.
  #'   \item \code{sd_cov}: reliability-weighted within-group standard
  #'     deviation.
  #'   \item \code{se_cov}: standard error of the mean.
  #'   \item \code{ci_low}, \code{ci_high}: lower/upper bounds of the
  #'     approximate 95\% confidence interval.
  #'   \item \code{n_samples}: number of curves contributing at this depth.
  #'   \item \code{n_observed}: curves contributing interpolated (observed)
  #'     data.
  #'   \item \code{n_projected}: curves contributing model-projected data.
  #' }
  #' Object attributes \code{group.levels}, \code{group.col}, \code{nps}, and
  #' \code{group} carry metadata for downstream use (e.g., by
  #' \code{plot.Nonpareil.Set} with \code{aggregate = TRUE}).
  #'
  #' @examples
  #' files <- system.file(
  #'   "extdata",
  #'   c("HumanGut.npo", "LakeLanier.npo", "IowaSoil.npo"),
  #'   package = "Nonpareil"
  #' )
  #' nps <- Nonpareil.set(files, plot = FALSE)
  #' # Two groups: first two curves vs. the third
  #' grp <- Nonpareil.group(nps, group = c("A", "A", "B"))
  #' head(grp)
  #'
  #' # Use plot.Nonpareil.Set convenience wrapper
  #' plot(nps, aggregate = TRUE, group = c("A", "A", "B"))
  invisible(out)
}

