#' Export Projection Surface as GeoTIFF
#'
#' Exports the projected illuminated surface as a GeoTIFF file with optional
#' spatial referencing.
#'
#' @param projection A projection_surface object from
#'   \code{\link{extract_projection_surface}}.
#' @param output_file Character string specifying the output file path.
#' @param crs Character or numeric. EPSG code or CRS string for spatial referencing
#'   (default: NULL for no spatial reference).
#' @param resolution Numeric. Pixel resolution in spatial units (default: 1).
#' @param origin Numeric vector of length 2 specifying the spatial origin
#'   (x, y coordinates). Default is c(0, 0).
#' @param compress Logical. Whether to use compression (default: TRUE).
#'
#' @details
#' The output GeoTIFF will contain:
#' - A single band with binary values (0 = not illuminated, 1 = illuminated)
#' - Spatial reference information (if crs is provided)
#' - Georeferencing information with specified resolution and origin
#'
#' @return Invisibly returns the path to the created file.
#'
#' @examples
#' \dontrun{
#'   export_projection_tif(
#'     projection = projection,
#'     output_file = "projection_output.tif",
#'     crs = "EPSG:4326"
#'   )
#' }
#'
#' @export
#' @importFrom terra rast writeRaster crs
export_projection_tif <- function(projection, output_file, crs = NULL,
                                 resolution = 1, origin = c(0, 0),
                                 compress = TRUE) {
  
  if (!inherits(projection, "projection_surface")) {
    stop("projection must be a projection_surface object")
  }
  
  if (!grepl("\\.tif$|\\.tiff$", output_file, ignore.case = TRUE)) {
    warning("Output file should have .tif or .tiff extension")
  }
  
  # Create parent directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get raster data
  raster_data <- projection$raster
  
  # Convert to terra raster
  r <- terra::rast(raster_data)
  
  # Set extent based on metadata and resolution
  x_range <- projection$metadata$x_range
  y_range <- projection$metadata$y_range
  
  # Create extent with proper origin
  ext <- terra::ext(
    origin[1] + x_range[1] * resolution,
    origin[1] + x_range[2] * resolution,
    origin[2] + y_range[1] * resolution,
    origin[2] + y_range[2] * resolution
  )
  
  terra::ext(r) <- ext
  
  # Set CRS if provided
  if (!is.null(crs)) {
    terra::crs(r) <- crs
  }
  
  # Write to file
  compression <- if (compress) "DEFLATE" else "NONE"
  
  tryCatch(
    {
      terra::writeRaster(r, output_file, overwrite = TRUE,
                        gdal = c("COMPRESS=" = compression))
      message("Projection exported to: ", output_file)
    },
    error = function(e) {
      stop("Failed to write GeoTIFF: ", e$message)
    }
  )
  
  invisible(output_file)
}


#' Export Reflectance Data as GeoTIFF (Multi-band)
#'
#' Exports selected spectral bands from the reflectance surface as a
#' multi-band GeoTIFF file.
#'
#' @param reflectance_data A reflectance_surface object.
#' @param output_file Character string specifying the output file path.
#' @param wavelengths Numeric vector of wavelengths to export.
#'   If NULL, exports all bands.
#' @param crs Character or numeric. EPSG code or CRS string.
#' @param resolution Numeric. Pixel resolution in spatial units.
#' @param origin Numeric vector of length 2 for spatial origin.
#'
#' @details
#' This function creates a multi-band GeoTIFF where each band corresponds to
#' a selected wavelength. Non-illuminated pixels are set to NA.
#'
#' @return Invisibly returns the path to the created file.
#'
#' @export
#' @importFrom terra rast writeRaster crs
export_reflectance_tif <- function(reflectance_data, output_file,
                                  wavelengths = NULL, crs = NULL,
                                  resolution = 1, origin = c(0, 0)) {
  
  if (!inherits(reflectance_data, "reflectance_surface")) {
    stop("reflectance_data must be a reflectance_surface object")
  }
  
  spectral_cube <- reflectance_data$spectral_cube
  illuminated_mask <- reflectance_data$illuminated_mask
  all_wavelengths <- reflectance_data$wavelengths
  
  # Select wavelengths
  if (is.null(wavelengths)) {
    band_indices <- seq_len(dim(spectral_cube)[3])
    selected_wavelengths <- all_wavelengths
  } else {
    band_indices <- numeric()
    selected_wavelengths <- numeric()
    
    for (wl in wavelengths) {
      idx <- which.min(abs(all_wavelengths - wl))
      band_indices <- c(band_indices, idx)
      selected_wavelengths <- c(selected_wavelengths, all_wavelengths[idx])
    }
  }
  
  # Create multi-band raster
  multi_band <- array(NA, dim = c(dim(spectral_cube)[1], dim(spectral_cube)[2],
                                 length(band_indices)))
  
  for (i in seq_along(band_indices)) {
    band_data <- spectral_cube[, , band_indices[i]]
    band_data[!illuminated_mask] <- NA
    multi_band[, , i] <- band_data
  }
  
  # Convert to terra raster
  r <- terra::rast(multi_band)
  
  # Set names
  names(r) <- paste0("band_", round(selected_wavelengths, 1), "nm")
  
  # Set extent
  proj <- reflectance_data$projection
  x_range <- proj$metadata$x_range
  y_range <- proj$metadata$y_range
  
  ext <- terra::ext(
    origin[1] + x_range[1] * resolution,
    origin[1] + x_range[2] * resolution,
    origin[2] + y_range[1] * resolution,
    origin[2] + y_range[2] * resolution
  )
  
  terra::ext(r) <- ext
  
  if (!is.null(crs)) {
    terra::crs(r) <- crs
  }
  
  # Write to file
  terra::writeRaster(r, output_file, overwrite = TRUE,
                    gdal = c("COMPRESS=" = "DEFLATE"))
  
  message("Reflectance data exported to: ", output_file)
  
  invisible(output_file)
}


#' Export Reflectance Statistics to CSV
#'
#' Saves mean reflectance statistics to a CSV file.
#'
#' @param reflectance_stats Data frame from \code{\link{get_mean_reflectance}}.
#' @param output_file Character string specifying the output CSV file path.
#' @param metadata List of additional metadata to include in the file header.
#'
#' @return Invisibly returns the path to the created file.
#'
#' @export
export_reflectance_csv <- function(reflectance_stats, output_file,
                                  metadata = NULL) {
  
  # Create parent directory if needed
  output_dir <- dirname(output_file)
  if (output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write metadata as comments if provided
  if (!is.null(metadata) && length(metadata) > 0) {
    meta_lines <- character()
    for (name in names(metadata)) {
      meta_lines <- c(meta_lines, paste0("# ", name, ": ", metadata[[name]]))
    }
    writeLines(meta_lines, output_file)
  }
  
  # Append the data
  write.csv(reflectance_stats, output_file, row.names = FALSE,
           append = !is.null(metadata))
  
  message("Reflectance statistics exported to: ", output_file)
  
  invisible(output_file)
}
