#' Single-cell RNA-seq Binarization Function
#'
#' This function binarize the single-cell celltag data
#' @param celltag.obj A CellTag object with the raw count matrix generated
#' @param tag.cutoff How many tags would you like to be used as a cutoff to say that the cells are tagged?
#' @return A CellTag object with the attribute (binary.mtx) filled.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' SingleCellDataBinatization(bam.test.obj, 2)
#' 
SingleCellDataBinatization <- function(celltag.obj, tag.cutoff) {
  CellTags <- celltag.obj@raw.count
  CellTags[CellTags < tag.cutoff] <- 0
  CellTags[CellTags > 0] <- 1
  celltag.obj@binary.mtx <- as(CellTags, "dgCMatrix")
  return(celltag.obj)
}

#' Single-cell RNA-seq Whitelisting Function
#'
#' This function conducts whitelist filtering through the single-cell dataset
#' @param celltag.obj A CellTag object with the binary matrix generated
#' @param whitels.cell.tag.file file director to the whitelisted cell tags
#' @return A CellTag object with the attribute (whitelisted.count) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' SingleCellDataWhitelist(bam.test.obj, "~/Desktop/My_Favourite_Whitelist.csv")
#' 
SingleCellDataWhitelist <- function(celltag.obj, whitels.cell.tag.file) {
  # Store the cell names
  CellTags <- celltag.obj@binary.mtx
  cell.names <- rownames(CellTags)
  
  # Process the celltag matrix to format below
  # row - celltag
  # col - cells
  CellTags <- t(CellTags)
  celltag.rownames <- row.names(CellTags)
  
  # Filter the matrix using whitelist
  if (endsWith(whitels.cell.tag.file, ".csv")) {
    separator <- ","
  } else {
    if (endsWith(whitels.cell.tag.file, ".txt") | endsWith(whitels.cell.tag.file, ".tsv")) {
      separator <- "\t"
    } else {
      separator <- " "
    }
  }
  whitelist <- read.delim(whitels.cell.tag.file, sep = separator, header = T, stringsAsFactors = F)
  whitelist.names <- whitelist[,1]
  whitelist <- Reduce(intersect, list(whitelist.names, celltag.rownames))
  celltags.whitelisted <- CellTags[whitelist,]
  colnames(celltags.whitelisted) <- cell.names
  
  celltag.obj@whitelisted.count <- as(as.matrix(celltags.whitelisted), "dgCMatrix")
  return(celltag.obj)
}


