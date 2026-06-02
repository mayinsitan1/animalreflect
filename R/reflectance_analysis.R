#' Extract Reflectance Data from Projection Surface
#'
#' Maps hyperspectral data to the projected surface of the 3D model.
#'
#' @param projection A projection_surface object from \code{\link{extract_projection_surface}}.
#' @param hyperspectral_data A hyperspectral_spe object from \code{\link{load_hyperspectral_spe}}.
#' @param wavelengths Numeric vector of wavelengths (in nm) corresponding to the spectral bands.
#' @param interpolation_method Character. Method for interpolating hyperspectral data
#'   to projection grid. Options: "nearest" (default), "bilinear".
#'
#' @details
#' This function aligns the hyperspectral cube with the projected surface and
#' extracts spectral information for each illuminated pixel. The output preserves
#' the full spectral dimension for later analysis.
#'
#' @return A list with class "reflectance_surface" containing:
#'   \item{spectral_cube}{3D array (x, y, bands) of reflectance values}
#'   \item{wavelengths}{Numeric vector of wavelengths (nm)}
#'   \item{illuminated_mask}{2D logical matrix of illuminated pixels}
#'   \item{projection}{Reference to the projection_surface object}
#'   \item{metadata}{List with processing information}
#'
#' @examples
#' \dontrun{
#'   reflectance <- extract_reflectance_surface(
#'     projection = projection,
#'     hyperspectral_data = hyper_data,
#'     wavelengths = hyper_data$wavelengths
#'   )
#' }
#'
#' @export
extract_reflectance_surface <- function(projection, hyperspectral_data,
                                       wavelengths, interpolation_method = "nearest") {
  
  if (!inherits(projection, "projection_surface")) {
    stop("projection must be a projection_surface object")
  }
  if (!inherits(hyperspectral_data, "hyperspectral_spe")) {
    stop("hyperspectral_data must be a hyperspectral_spe object")
  }
  
  if (is.null(hyperspectral_data$data)) {
    stop("Hyperspectral data is empty. Check if SPE file was properly loaded.")
  }
  
  if (!interpolation_method %in% c("nearest", "bilinear")) {
    stop("interpolation_method must be 'nearest' or 'bilinear'")
  }
  
  # Get dimensions
  proj_raster <- projection$raster
  proj_dim <- dim(proj_raster)
  hyper_dim <- dim(hyperspectral_data$data)
  
  # Create spectral cube for projected surface
  spectral_cube <- array(NA, dim = c(proj_dim[1], proj_dim[2], length(wavelengths)))
  
  # Map hyperspectral data to projection grid
  if (interpolation_method == "nearest") {
    # Nearest neighbor interpolation
    for (i in seq_len(proj_dim[1])) {
      for (j in seq_len(proj_dim[2])) {
        if (proj_raster[i, j] == 0) next  # Skip non-illuminated pixels
        
        # Map projection coordinates to hyperspectral image coordinates
        hyper_i <- round((i - 1) / proj_dim[1] * hyper_dim[1]) + 1
        hyper_j <- round((j - 1) / proj_dim[2] * hyper_dim[2]) + 1
        
        hyper_i <- max(1, min(hyper_i, hyper_dim[1]))
        hyper_j <- max(1, min(hyper_j, hyper_dim[2]))
        
        spectral_cube[i, j, ] <- hyperspectral_data$data[hyper_i, hyper_j, ]
      }
    }
  } else if (interpolation_method == "bilinear") {
    # Bilinear interpolation (simplified)
    for (i in seq_len(proj_dim[1])) {
      for (j in seq_len(proj_dim[2])) {
        if (proj_raster[i, j] == 0) next
        
        # Calculate fractional coordinates
        hyper_i_frac <- (i - 1) / proj_dim[1] * hyper_dim[1] + 1
        hyper_j_frac <- (j - 1) / proj_dim[2] * hyper_dim[2] + 1
        
        hyper_i1 <- floor(hyper_i_frac)
        hyper_i2 <- ceiling(hyper_i_frac)
        hyper_j1 <- floor(hyper_j_frac)
        hyper_j2 <- ceiling(hyper_j_frac)
        
        hyper_i1 <- max(1, min(hyper_i1, hyper_dim[1]))
        hyper_i2 <- max(1, min(hyper_i2, hyper_dim[1]))
        hyper_j1 <- max(1, min(hyper_j1, hyper_dim[2]))
        hyper_j2 <- max(1, min(hyper_j2, hyper_dim[2]))
        
        # Interpolation weights
        wi <- hyper_i_frac - hyper_i1
        wj <- hyper_j_frac - hyper_j1
        
        # Bilinear interpolation
        val1 <- (1 - wi) * (1 - wj) * hyperspectral_data$data[hyper_i1, hyper_j1, ]
        val2 <- wi * (1 - wj) * hyperspectral_data$data[hyper_i2, hyper_j1, ]
        val3 <- (1 - wi) * wj * hyperspectral_data$data[hyper_i1, hyper_j2, ]
        val4 <- wi * wj * hyperspectral_data$data[hyper_i2, hyper_j2, ]
        
        spectral_cube[i, j, ] <- val1 + val2 + val3 + val4
      }
    }
  }
  
  # Create illuminated mask
  illuminated_mask <- proj_raster > 0
  
  # Create result object
  reflectance_data <- list(
    spectral_cube = spectral_cube,
    wavelengths = wavelengths,
    illuminated_mask = illuminated_mask,
    projection = projection,
    metadata = list(
      n_illuminated_pixels = sum(illuminated_mask),
      n_bands = length(wavelengths),
      wavelength_range = c(min(wavelengths), max(wavelengths)),
      interpolation_method = interpolation_method
    )
  )
  
  class(reflectance_data) <- "reflectance_surface"
  
  return(reflectance_data)
}


#' Extract Mean Reflectance at Specific Wavelengths
#'
#' Calculates mean reflectance values at specified wavelengths from the
#' extracted reflectance surface.
#'
#' @param reflectance_data A reflectance_surface object from
#'   \code{\link{extract_reflectance_surface}}.
#' @param wavelengths Numeric vector of wavelengths (in nm) to extract.
#'   Must be within the range of available wavelengths.
#' @param hyperspectral_wavelengths Numeric vector of all available wavelengths
#'   in the hyperspectral data.
#' @param method Character. Aggregation method: "mean" (default), "median", or "all".
#'   "all" returns values for all illuminated pixels.
#'
#' @details
#' This function interpolates the spectral cube to the specified wavelengths
#' (if they don't exactly match available bands) and computes statistics.
#'
#' @return A data frame (or list if method="all") containing:
#'   \item{wavelength}{The specified wavelengths (nm)}
#'   \item{mean_reflectance}{Mean reflectance value across illuminated pixels}
#'   \item{std_dev}{Standard deviation of reflectance}
#'   \item{min}{Minimum reflectance value}
#'   \item{max}{Maximum reflectance value}
#'
#' @examples
#' \dontrun{
#'   # Get mean reflectance at specific wavelengths
#'   mean_refl <- get_mean_reflectance(
#'     reflectance_data = reflectance,
#'     wavelengths = c(550, 850, 1650)
#'   )
#'   print(mean_refl)
#' }
#'
#' @export
get_mean_reflectance <- function(reflectance_data, wavelengths,
                                hyperspectral_wavelengths = NULL,
                                method = "mean") {
  
  if (!inherits(reflectance_data, "reflectance_surface")) {
    stop("reflectance_data must be a reflectance_surface object")
  }
  
  if (is.null(hyperspectral_wavelengths)) {
    hyperspectral_wavelengths <- reflectance_data$wavelengths
  }
  
  if (!method %in% c("mean", "median", "all")) {
    stop("method must be 'mean', 'median', or 'all'")
  }
  
  spectral_cube <- reflectance_data$spectral_cube
  illuminated_mask <- reflectance_data$illuminated_mask
  
  results <- data.frame(
    wavelength = numeric(),
    mean_reflectance = numeric(),
    std_dev = numeric(),
    min = numeric(),
    max = numeric()
  )
  
  for (wl in wavelengths) {
    # Find closest bands in the spectral cube
    idx <- which.min(abs(hyperspectral_wavelengths - wl))
    
    if (abs(hyperspectral_wavelengths[idx] - wl) > 5) {
      warning("Wavelength ", wl, " not closely matched. Closest is ",
              hyperspectral_wavelengths[idx])
    }
    
    # Extract reflectance values for this wavelength
    band_data <- spectral_cube[, , idx]
    illuminated_values <- band_data[illuminated_mask]
    
    # Remove NA values
    illuminated_values <- illuminated_values[!is.na(illuminated_values)]
    
    if (length(illuminated_values) == 0) {
      warning("No valid data for wavelength ", wl)
      next
    }
    
    if (method == "mean") {
      results <- rbind(results, data.frame(
        wavelength = wl,
        mean_reflectance = mean(illuminated_values, na.rm = TRUE),
        std_dev = sd(illuminated_values, na.rm = TRUE),
        min = min(illuminated_values, na.rm = TRUE),
        max = max(illuminated_values, na.rm = TRUE)
      ))
    } else if (method == "median") {
      results <- rbind(results, data.frame(
        wavelength = wl,
        mean_reflectance = median(illuminated_values, na.rm = TRUE),
        std_dev = sd(illuminated_values, na.rm = TRUE),
        min = min(illuminated_values, na.rm = TRUE),
        max = max(illuminated_values, na.rm = TRUE)
      ))
    }
  }
  
  if (method == "all") {
    return(list(
      wavelengths = wavelengths,
      values = illuminated_values
    ))
  }
  
  return(results)
}


#' Print method for reflectance_surface
#'
#' @param x A reflectance_surface object
#' @param ... Additional arguments passed to print
#'
#' @export
print.reflectance_surface <- function(x, ...) {
  cat("Reflectance Surface (Hyperspectral)\n")
  cat("===================================\n")
  cat("Illuminated pixels:", x$metadata$n_illuminated_pixels, "\n")
  cat("Spectral bands:", x$metadata$n_bands, "\n")
  cat("Wavelength range:", round(x$metadata$wavelength_range[1], 1), "-",
      round(x$metadata$wavelength_range[2], 1), "nm\n")
  cat("Interpolation method:", x$metadata$interpolation_method, "\n")
}
