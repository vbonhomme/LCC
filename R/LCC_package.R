##### Package documentation and NAMESPACE import

#' LCC (Large Charcoal Count)
#'
#' Here comes the Description
#'
#' @docType package
#' @name LCC_package
#' @import dplyr
NULL

# prevents "no visible binding for global variable"
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
globalVariables(c(
  "Seedle",
  "Age_calBP",
  "Area",
  "Count",
  "cpeak",
  "FRI",
  "pic",
  "Pic",
  "Sarea",
  "SdlArea",
  "SdlCounts",
  "Se.area.inf",
  "Se.area.sup",
  "Sinf",
  "Smpl",
  "Ssup",
  "Suminf",
  "Sumsup",
  "thresh"))
