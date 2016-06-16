library(WGCNA)
library(flashClust)
library(matrixcalc)
#####################################
heatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                       distfun = dist, hclustfun = hclust, dendrogram = c("both",
                                                                          "row", "column", "none"), symm = FALSE, scale = c("none",
                                                                                                                            "row", "column"), na.rm = TRUE, revC = identical(Colv,
                                                                                                                                                                             "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
                         scale != "none", col = "heat.colors", colsep, rowsep,
                       sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
                       notecol = "cyan", na.color = par("bg"), trace = c("column",
                                                                         "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
                       vline = median(breaks), linecol = tracecol, margins = c(5,
                                                                               5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
                       cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
                       key = TRUE, keysize = 1.5, density.info = c("histogram",
                                                                   "density", "none"), denscol = tracecol, symkey = min(x <
                                                                                                                          0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
                       xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, ...)
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) !=
            nc)
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                      1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) !=
            nr)
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                                           1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1],
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



heatmap.3 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                       distfun = dist, hclustfun = hclust, dendrogram = c("both",
                                                                          "row", "column", "none"), symm = FALSE, scale = c("none",
                                                                                                                            "row", "column"), na.rm = TRUE, revC = identical(Colv,
                                                                                                                                                                             "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
                         scale != "none", col = "heat.colors", colsep, rowsep,
                       sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
                       notecol = "cyan", na.color = par("bg"), trace = c("column",
                                                                         "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
                       vline = median(breaks), linecol = tracecol, margins = c(5,
                                                                               5), ColSideColors, RowSideColors, cexRow = 2*(0.2 + 1/log10(nr)),
                       cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
                       key = TRUE, keysize = 1.5, density.info = c("histogram",
                                                                   "density", "none"), denscol = tracecol, symkey = min(x <
                                                                                                                          0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
                       xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
                       ...)
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }

  if (exists("hcc")) {retval$hcc <- hcc}
  if (exists("hcr")) {retval$hcr <- hcr}
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors)) #|| ncol(ColSideColors) != nc)
        stop("'ColSideColors' must be a character ") #vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      if (is.vector(ColSideColors)) nnn=1
      else nnn=nrow(ColSideColors)
      lhei <- c(lhei[1], nnn*0.1, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors)) #|| length(RowSideColors) != nr)
        stop("'RowSideColors' must be a character ")
      if (is.vector(RowSideColors)) nnn=1
      else nnn=ncol(RowSideColors)
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], nnn*0.1, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))

    if (is.vector(RowSideColors)) {
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (is.matrix(RowSideColors)) {
      jk.row = RowSideColors
      jk.xy = matrix(which(jk.row != "0"), dim(jk.row))
      colnames(jk.xy) <- colnames(jk.row)
      #image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
      #grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
      image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
      #            axis(3, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
      axis(1, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
      #   grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
    }

  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    if (is.vector(ColSideColors)) {
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    if (is.matrix(ColSideColors)) {
      jk.col = ColSideColors
      jk.xy = matrix(which(jk.col != "0"), dim(jk.col))
      image(t(jk.xy), col = jk.col[, colInd], axes = FALSE)
    }

  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
    cellnote <- t(cellnote)
  }
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  # if (!invalid(na.color) & any(is.na(x))) {
  #     mmat <- ifelse(is.na(x), 1, NA)
  #    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
  #         col = na.color, add = TRUE)
  # }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  #    axis(3, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
  #        cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 3, line = margins[1] - 1.25)

  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1],
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x))
      tmpbreaks[length(tmpbreaks)] <- max(abs(x))
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else {
      #mtext(side = 1, "Value", line = 2)
      mtext(side = 1, "", line = 2)
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("") ;#Color Key and Density Plot", cex=0.25)
      #title("Color Key and Density Plot", cex=0.25)
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
      #mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      #title("Color Key and Histogram", cex=0.25)
      title("", cex=0.25)
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
      #mtext(side = 2, "Count", line = 2)
    }
    else { title("", cex=0.25)
           #title("Color Key", cex=0.25)
    }
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

test.cv.rf <- function(datver="09102015", norm.dat, cl, n.bin=5,padj.th = 0.01, n.markers=10,othervals=c()) {

  require(randomForest)
  require(e1071)
  bins = sample(1:n.bin, length(cl),replace=T)
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    print(c(length(x), length(tmp)))
    setNames(tmp[sample(length(tmp))], x)
  }))

  if (datver=="09102015") {
    tmp1 <- get.field(names(bins), ";", 1)
    tmp12 <- get.field(tmp1, ".", 2)
    tmp2 <- get.field(names(bins), ";", 2)
    names(bins) <- paste(tmp12, tmp2, sep=";")
  } else {
    tmp1 <- get.field(names(bins), ";", 1)
    tmp2 <- get.field(names(bins), ";", 2)
    tmp22 <- get.field(tmp2, "-", 2)
    names(bins) <- paste(tmp1, tmp22, sep=";")
  }
  #names(cl) <- paste(cl,names(cl), sep=".")
  #bins <- bins[names(cl)]
  #print(bins)
  rf.pred.cl = setNames(rep(NA, length(cl)), names(cl))
  rf.pred.prob = matrix(0, nrow=length(cl), ncol=length(unique(cl)))
  rownames(rf.pred.prob) <- names(cl)
  colnames(rf.pred.prob) <- unique(cl)
  elecmp2<-c()
  for(i in 1:n.bin){
    if (n.bin==1) {
      select.cells = names(cl)
    } else {
      select.cells = names(cl)[bins!=i]
    }


    de.df = DE.genes.pw(norm.dat[,select.cells], paste("cl",cl[select.cells],sep=""))
    de.genes = sapply(de.df, function(x){ x = x[order(x$padj),]
                                          up = x$padj < padj.th & x$lfc > 1
                                          down = x$padj < padj.th & x$lfc < -1
                                          c(head(row.names(x)[up],n.markers), head(row.names(x)[down],n.markers))
    }, simplify=F)
    #markers = as.numeric(unique(unlist(de.genes)))
    markers = unique(unlist(de.genes))
    if (length(markers)>1) {
      rf.result<-randomForest(as.matrix(t(norm.dat[markers, select.cells])),as.factor(cl[select.cells]),ntree=1000)

      tmp = predict(rf.result, as.matrix(t(norm.dat[markers, names(cl)[bins==i]])),type="prob")
      rf.pred.prob[bins==i,colnames(tmp)] <- tmp
      rf.pred.cl[bins ==i] <-  apply(tmp, 1, which.max)
      if (length(othervals)>0) {
        tmp2 = predict(rf.result, as.matrix(t(norm.dat[markers, names(othervals)])),type="prob")
        tmp2 = apply(tmp2,1,which.max)
      }
    } else {
      rf.pred.cl[bins == i] = NA;
    }
  }
  if ((length(othervals)>0) & (length(tmp2)>0)) {
    names(tmp2)<-names(othervals)
  } else {
    tmp2<-rep(NA,length(othervals))
    names(tmp2)<-names(othervals)
  }
  return(list(binvals=bins,rf=list(class=rf.pred.cl, posterior = rf.pred.prob,classother=tmp2)))
}

DE.genes.pw <- function(dat,cl){
  require(limma)
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  de.df=list()
  for(i in 1:(length(tmp.cl)-1)){
    for(j in (i+1):length(tmp.cl)){
      x=tmp.cl[i]
      y=tmp.cl[j]
      ctr <<- paste(x, "- ", y)
      contrasts.matrix <- makeContrasts(ctr,  levels=design)
      fit2 = contrasts.fit(fit, contrasts.matrix)
      fit2 = eBayes(fit2)
      padj = apply(fit2$p.value, 2, p.adjust)
      lfc = coef(fit2)
      pair=paste(x,y,sep="_")
      de.df[[pair]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1])
      row.names(de.df[[pair]])= row.names(dat)
    }
  }
  return(de.df)
}

###Function to run cross-validation on clusters using Random Forest. This code divides
###the overall sample set into groups of 20%, and then fits a random forest to the 80%
###and predicts the membership of the remaining 20%. All comparisons are pairwise among
###clusters (to minimize limma's inclusion of non-discriminatory genes in DE analysis)
###and then non-dominated clusters are assigned to each sample.
RUN_CVRF <- function ( AssignedID_FN, Feature, outdir, flag.plot=TRUE, datver="09102015", numruns=100 ) {

  #numruns<-100;  ###number of cross-validation runs
  pthresh<-0.05  ###threshold p-value for differentially expressed genes used for classification
  keepmarkers<-10  ###number of differentially expressed genes used for classification among each pair
  outfile <- paste("RFsummary.full_assignment_rf_mc.", numruns, ".10fold.RData", sep="")   ###output Rdata file with all cross-validation info
  outfile2<- paste("RFsummary.classification_primary_secondary_2_mc.", numruns, ".10fold.csv", sep="")  ###output csv with primary and secondary membership
  outfile3<- paste("RFsummary.10run_mc.", numruns, ".10fold.pdf", sep="")  ###output csv with primary and secondary membership
  outfile4<- paste("RFsummary.classification_Any.Same.Assignment.", numruns, ".10fold.csv", sep="")  ###output csv with primary and secondary membership


  filename.Id <- AssignedID_FN
  out <- read.csv(filename.Id, header=TRUE)
  Id.name <- gsub("~sim", "",as.character(out[,1]))
  Id.cre <- get.field(as.character(out[,1]), ";",1)
  Id.id <- as.character(out[,2])
  names(Id.id) <- Id.name

  XXX <- Feature
  #rownames(XXX) <- gsub("-197628.06.01.01;NA", ";197628.06.01.01", gsub("-197628.06.02.01;NA", ";197628.06.02.01", gsub("-197628.03.02.01;NA", ";197628.03.02.01", rownames(XXX))))
  nfeat <- ncol(XXX)
  idx <- match(rownames(XXX), names(Id.id))

  XXX.ClusterID <- Id.id[idx]
  unique.ClusterID <- sort(unique(XXX.ClusterID))
  XXX.label <- match(XXX.ClusterID, unique.ClusterID)
  names(XXX.label) <- names(Id.id)[idx]
  unique.label <- unique(XXX.label)

  set.seed(0) ###random number seed
  #excludeclusters<-c(0,11,13,2,4,7)  ###clusters to exclude in the cross-validation, first round
  excludeclusters<-c() ###clusters to exclude in the cross-validation, second round
  cross.validation=T;  ###T=cross-validation (i.e. removing 20% of cells), F=classification of new data

  rpkmcount<-t(XXX)
  allcl<-paste0("c",XXX.label)
  names(allcl) <- names(XXX.label)

  shuffmat<-c();
  for (ii in (unique(allcl))) {
    shuffmat<-cbind(shuffmat,rpkmcount[,names(allcl)[which(allcl==ii)[1]]])
  }
  rpkmcount2<-cbind(rpkmcount,shuffmat)
  tempvec2<-rep("Out",ncol(shuffmat))
  names(tempvec2)<-paste("Out",1:ncol(shuffmat),sep=";")
  allcl<-c(allcl,tempvec2)
  clusttab<-table(allcl)
  allclust<-setdiff(unique(allcl),paste0("c",excludeclusters))
  allclust<-allclust[order(allclust)]
  allcells<-names(allcl[allcl!="Out"])
  colnames(rpkmcount2)<-names(allcl)

  ####remove 20% of cells, run cross-validation against every pair of clusters###
  full_assignment_rf_all<-list();
  if (cross.validation) {
    nbins<-10;
  } else {
    nbins<-1;
  }

  for (ii in 1:numruns) {
    full_assignment_rf_all[[ii]] <-test.cv.rf(datver, rpkmcount2[,],allcl,n.bin=nbins,padj.th=pthresh,n.markers=keepmarkers,othervals=c())
  }


  for (ii in 1:numruns) {
    if (ii==1) {
      Ncluster <- length(unique(allcl))
      NN <- length(full_assignment_rf_all[[ii]]$rf$class)
      fulltab <- matrix(0, nrow=NN, ncol=Ncluster)
      rownames(fulltab) <- names(full_assignment_rf_all[[ii]]$rf$class)
      colnames(fulltab) <- sort(unique(allcl))
    }

    ii.names <- names(full_assignment_rf_all[[ii]]$rf$class)
    ii.class <- full_assignment_rf_all[[ii]]$rf$class
    for (n in 1:NN) {
      fulltab[ii.names[n], ii.class[n]] <- fulltab[ii.names[n], ii.class[n]]  + 1
    }
  }

  colnames(fulltab) <- paste0("c",unique.ClusterID[as.numeric(gsub("c","",colnames(fulltab)))])
  nc <- length(colnames(fulltab))
  nx <- length(XXX.ClusterID)
  reordered <- c(order(colnames(fulltab)[1:(nc-1)]) , nc)
  fulltab.reordered <- fulltab[1:nx, reordered]


  save(full_assignment_rf_all,fulltab, fulltab.reordered, file=paste(outdir, outfile, sep="/"))


  fulltab.reordered[fulltab.reordered > 100] <- 100
  primvec <- apply(fulltab.reordered, 1, function(x) { order(x,decreasing=TRUE)[c(1)] })
  secondvec <- apply(fulltab.reordered, 1, function(x) { order(x,decreasing=TRUE)[c(2)] })
  outdat<-data.frame(cell=allcells,primvec=primvec,secondvec=secondvec,fulltab.reordered)
  outdat<-outdat[outdat$cell %in% colnames(rpkmcount2),]
  write.csv(outdat,file=paste(outdir, outfile2, sep="/"))

  creline_table = read.csv(paste0(data_dir,"Cre_line_typeC.09102015.csv"), header=TRUE)
  creC <- as.character(creline_table[,"cre_line"])
  mycre <- get.field(rownames(fulltab.reordered), ";",1)
  crecolor <- rainbow(11+1)[match(mycre,creC)]


  leg.pch <- rep(15,11)
  leg.str <- creC
  leg.crecolor <- rainbow(12)[1:11]

  myheat.color <- rev(heat.colors(101)) ; myheat.color[101] <- "black"
  if (flag.plot) {
    pdf(paste(outdir, outfile3, sep="/"), height=8, width=16)
    heatmap.3(t(fulltab.reordered), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main="10-fold cross validation (100 runs) ")
    legend("bottomleft", col=leg.crecolor, pch=rep(15,11), legend=leg.str,cex=0.75)
  }
  print(" summary table")
  print(" fulltab.reordered  (702 x Ncluster) ")
  idx.cNA <- match("cNA", colnames(fulltab.reordered))
  print("Samples Always Assigned to NA")
  print(rownames(fulltab.reordered)[fulltab.reordered[, idx.cNA]==10])

  fulltab.reordered.woNA <- fulltab.reordered[, -idx.cNA]

  tmpR <- matrix(rep(apply(fulltab.reordered.woNA, 1, sum), ncol(fulltab.reordered.woNA)), nrow=nrow(fulltab.reordered.woNA))
  fulltab.reordered.woNA.R <- round(100 * fulltab.reordered.woNA / tmpR)
  fulltab.reordered.woNA.R[ is.na(fulltab.reordered.woNA.R) ] <- 0
  class.N2R <- t(apply(fulltab.reordered.woNA.R, 2, function(x) { Ngt0 <- sum(x>0); N100 <- sum(x==100) ; c(Ngt0,N100,round(100*N100/Ngt0)) } ))
  colnames(class.N2R) <- c("N.Any", "N.Same", "Same(%)")
  print(class.N2R)
  write.csv(class.N2R, file=paste(outdir, outfile4, sep="/"))

  #colnames(fulltab.reordered.woNA.R) <- paste(class.N2R[,3],colnames(fulltab.reordered.woNA.R),sep="%")
  colnames(fulltab.reordered.woNA.R) <- paste(colnames(fulltab.reordered.woNA.R),' (',class.N2R[,3],'%)')
  if (flag.plot) {
    heatmap.3(t(fulltab.reordered.woNA.R), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main="100 run of 10-fold cross validation")
    legend("bottomleft", col=leg.crecolor, pch=rep(15,11), legend=leg.str,cex=0.75)
    dev.off()
  }
}















