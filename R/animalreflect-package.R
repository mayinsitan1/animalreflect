#' animalreflect: Integration of 3D Animal Models with Hyperspectral Imaging
#'
#' An R package for integrating 3D animal models with hyperspectral imaging data
#' to extract projection surfaces based on solar geometry and retrieve spectral
#' reflectance information.
#'
#' @details
#' The package provides tools for:
#' \itemize{
#'   \item Loading 3D animal models in OBJ format
#'   \item Loading hyperspectral images in SPE format
#'   \item Calculating solar position from date, time, and location
#'   \item Extracting illuminated projection surfaces
#'   \item Retrieving spectral reflectance data
#'   \item Exporting results as GeoTIFF and CSV files
#'   \item Batch processing multiple models
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{load_obj_model}}: Load 3D OBJ model files
#'   \item \code{\link{load_hyperspectral_spe}}: Load SPE format hyperspectral images
#'   \item \code{\link{calculate_solar_position}}: Calculate solar position
#'   \item \code{\link{extract_projection_surface}}: Extract illuminated projection
#'   \item \code{\link{extract_reflectance_surface}}: Extract reflectance data
#'   \item \code{\link{get_mean_reflectance}}: Calculate mean reflectance statistics
#'   \item \code{\link{export_projection_tif}}: Export projection as GeoTIFF
#'   \item \code{\link{batch_process_animals}}: Batch process multiple models
#' }
#'
#' @section Data Formats:
#' \describe{
#'   \item{OBJ}{Wavefront 3D model format}
#'   \item{SPE}{Princeton Instruments hyperspectral image format}
#'   \item{GeoTIFF}{Georeferenced raster output format}
#'   \item{CSV}{Reflectance statistics export format}
#' }
#'
#' @section Getting Started:
#' See the README and vignettes for detailed examples and workflows.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom terra rast writeRaster crs ext
#' @importFrom lubridate ymd_hms as_datetime with_tz hour minute second
#'
"_PACKAGE"

## usethis namespace: start
#' @useDynLib animalreflect, .registration = TRUE
## usethis namespace: end
NULL
