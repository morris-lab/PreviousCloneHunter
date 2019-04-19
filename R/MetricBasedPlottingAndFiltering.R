#' Metric-Base Filtering Function
#'
#' This function applies further filtering on scRNA-seq data with CellTags based on cutoff values identified from the metric plots
#' @param celltag.obj A CellTag Object with count matrix generated
#' @param cutoff The cutoff decided from the metric plots
#' @param comparison Would you like to maintain the part less than/greater than the cutoff? Default to less. Choices can be greater or less.
#' @return A CellTag Object with attribute (metric.filtered.count) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricBasedFiltering(bam.test.object, 20, "less")
#'
MetricBasedFiltering <- function(celltag.obj, cutoff, comparison = "less") {
  whitelisted.celltag.data <- t(celltag.obj@whitelisted.count)
  # Set up the filtering data frame
  CellTags.per.cell.whitelisted.pf <- as.data.frame(rowSums(whitelisted.celltag.data))
  
  # Set up the filtered celltag dataset object
  if (comparison == "less") {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf < (cutoff + 1))
  } else {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf > (cutoff - 1))
  }
  cell.bc.filter <- row.names(cell.filter)
  # Filter celltag dataset
  celltags.whitelisted.new <- whitelisted.celltag.data[cell.bc.filter, ]
  
  celltag.obj@metric.filtered.count <- celltags.whitelisted.new
  return(celltag.obj)
}

#' CellTag Metric Plotting Function
#'
#' This function provides some metric plots for further downstream celltag filtering in the scRNA-seq dataset
#' @param celltag.obj A CellTag Object
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricPlots(bam.test.obj)
#'
MetricPlots <- function(celltag.obj) {
  if (nrow(celltag.obj@metric.filtered.count) <= 0) {
    if (nrow(celltag.obj@whitelisted.count) <= 0) {
      celltag.data <- celltag.obj@binary.mtx
    } else {
      celltag.data <- celltag.obj@whitelisted.count
    }
  } else {
    celltag.data <- celltag.obj@metric.filtered.count
  }
  
  CellTags.per.cell.whitelisted.pf <- rowSums(celltag.data)
  CellTags.per.cell.avg <- mean(CellTags.per.cell.whitelisted.pf)
  CellTags.frequency.whitelisted.pf <- colSums(celltag.data)
  CellTags.freq.avg <- mean(CellTags.frequency.whitelisted.pf)
  plot(CellTags.per.cell.whitelisted.pf)
  plot(CellTags.frequency.whitelisted.pf)
  cat("Average: ", CellTags.per.cell.avg, "\n")
  cat("Frequency: ", CellTags.freq.avg, "\n")
}

