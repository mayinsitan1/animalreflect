#' Batch Process Multiple Animal Models
#'
#' Processes multiple 3D models and hyperspectral images in batch mode with
#' consistent solar geometry and export settings.
#'
#' @param model_files Character vector of paths to OBJ model files.
#' @param hyperspectral_files Character vector of paths to SPE hyperspectral files.
#' @param date Character string in format "YYYY-MM-DD".
#' @param time Character string in format "HH:MM:SS".
#' @param latitude Numeric. Latitude for solar calculation.
#' @param longitude Numeric. Longitude for solar calculation.
#' @param timezone Character. Time zone for calculations.
#' @param output_dir Character. Directory for output files.
#' @param export_projection Logical. Whether to export projection surfaces (default: TRUE).
#' @param export_reflectance Logical. Whether to export reflectance data (default: FALSE).
#' @param wavelengths_extract Numeric vector of wavelengths for reflectance extraction.
#' @param crs Character. EPSG code or CRS string for output files.
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @details
#' This function performs the following for each model:
#' 1. Loads OBJ model
#' 2. Loads hyperspectral image
#' 3. Calculates solar position
#' 4. Extracts projection surface
#' 5. Extracts reflectance data (if enabled)
#' 6. Exports results to specified formats
#'
#' Results are named based on the input filenames and stored in output_dir.
#'
#' @return A data frame with processing results for each model containing:
#'   \item{model_file}{Input model filename}
#'   \item{hyper_file}{Input hyperspectral filename}
#'   \item{status}{Processing status ("success" or error message)}
#'   \item{projection_file}{Output projection TIF path (if exported)}
#'   \item{reflectance_file}{Output reflectance TIF path (if exported)}
#'   \item{n_illuminated_pixels}{Number of illuminated pixels}
#'   \item{processing_time}{Time taken for processing (seconds)}
#'
#' @examples
#' \dontrun{
#'   models <- list.files("data/models", pattern = "\\.obj$", full.names = TRUE)
#'   hyperspectral <- list.files("data/images", pattern = "\\.spe$", full.names = TRUE)
#'
#'   results <- batch_process_animals(
#'     model_files = models,
#'     hyperspectral_files = hyperspectral,
#'     date = "2026-06-01",
#'     time = "12:00:00",
#'     latitude = 40.7128,
#'     longitude = -74.0060,
#'     output_dir = "output/results",
#'     wavelengths_extract = c(550, 850, 1650)
#'   )
#'
#'   print(results)
#' }
#'
#' @export
batch_process_animals <- function(model_files, hyperspectral_files,
                                 date, time, latitude, longitude,
                                 timezone = "UTC", output_dir = "./output",
                                 export_projection = TRUE,
                                 export_reflectance = FALSE,
                                 wavelengths_extract = NULL,
                                 crs = NULL, verbose = TRUE) {
  
  # Validation
  if (length(model_files) != length(hyperspectral_files)) {
    stop("model_files and hyperspectral_files must have the same length")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  n_files <- length(model_files)
  results <- data.frame(
    model_file = character(n_files),
    hyper_file = character(n_files),
    status = character(n_files),
    projection_file = character(n_files),
    reflectance_file = character(n_files),
    n_illuminated_pixels = numeric(n_files),
    processing_time = numeric(n_files),
    stringsAsFactors = FALSE
  )
  
  # Calculate solar position once (same for all)
  if (verbose) message("Calculating solar position for ", date, " ", time)
  
  solar_position <- calculate_solar_position(
    date = date,
    time = time,
    latitude = latitude,
    longitude = longitude,
    timezone = timezone
  )
  
  if (verbose) {
    cat("Solar geometry: Azimuth =", round(solar_position$azimuth, 1), "°,",
        "Elevation =", round(solar_position$elevation, 1), "°\n")
  }
  
  # Process each file pair
  for (i in seq_len(n_files)) {
    start_time <- Sys.time()
    
    model_file <- model_files[i]
    hyper_file <- hyperspectral_files[i]
    
    results$model_file[i] <- basename(model_file)
    results$hyper_file[i] <- basename(hyper_file)
    
    if (verbose) {
      cat("\n[", i, "/", n_files, "] Processing ",
          basename(model_file), " with ", basename(hyper_file), "...\n", sep = "")
    }
    
    tryCatch({
      # Load model
      if (verbose) cat("  Loading model...\n")
      model <- load_obj_model(model_file)
      
      # Load hyperspectral data
      if (verbose) cat("  Loading hyperspectral data...\n")
      hyper_data <- load_hyperspectral_spe(hyper_file)
      
      # Extract projection
      if (verbose) cat("  Extracting projection surface...\n")
      projection <- extract_projection_surface(model, solar_position)
      
      results$n_illuminated_pixels[i] <- projection$metadata$n_illuminated_faces
      
      # Export projection if requested
      if (export_projection) {
        base_name <- tools::file_path_sans_ext(basename(model_file))
        proj_file <- file.path(output_dir, paste0(base_name, "_projection.tif"))
        
        if (verbose) cat("  Exporting projection to", basename(proj_file), "\n")
        export_projection_tif(projection, proj_file, crs = crs)
        results$projection_file[i] <- proj_file
      }
      
      # Extract and export reflectance if requested
      if (export_reflectance && !is.null(hyper_data$data)) {
        if (verbose) cat("  Extracting reflectance data...\n")
        reflectance <- extract_reflectance_surface(
          projection, hyper_data,
          hyper_data$wavelengths
        )
        
        base_name <- tools::file_path_sans_ext(basename(model_file))
        refl_file <- file.path(output_dir, paste0(base_name, "_reflectance.tif"))
        
        if (verbose) cat("  Exporting reflectance to", basename(refl_file), "\n")
        export_reflectance_tif(
          reflectance, refl_file,
          wavelengths = wavelengths_extract,
          crs = crs
        )
        results$reflectance_file[i] <- refl_file
      }
      
      results$status[i] <- "success"
      
      if (verbose) cat("  ✓ Completed successfully\n")
      
    }, error = function(e) {
      results$status[i] <<- paste("ERROR:", e$message)
      if (verbose) cat("  ✗ Error:", e$message, "\n")
    })
    
    results$processing_time[i] <- as.numeric(Sys.time() - start_time, units = "secs")
  }
  
  # Summary
  if (verbose) {
    cat("\n========== BATCH PROCESSING SUMMARY ==========\n")
    cat("Total files processed:", n_files, "\n")
    cat("Successful:", sum(results$status == "success"), "\n")
    cat("Failed:", sum(results$status != "success"), "\n")
    cat("Total processing time:", sum(results$processing_time), "seconds\n")
  }
  
  return(results)
}


#' Create Processing Report
#'
#' Generates a summary report of batch processing results.
#'
#' @param batch_results Data frame from \code{\link{batch_process_animals}}.
#' @param output_file Optional path to save the report as text file.
#'
#' @return A formatted character vector (invisible) summarizing the results.
#'
#' @export
create_processing_report <- function(batch_results, output_file = NULL) {
  
  report_lines <- c(
    "ANIMAL REFLECTANCE ANALYSIS - PROCESSING REPORT",
    "=============================================\n",
    paste("Processing Date:", Sys.time()),
    paste("Total Models Processed:", nrow(batch_results)),
    paste("Successful:", sum(batch_results$status == "success")),
    paste("Failed:", sum(batch_results$status != "success")),
    paste("Total Processing Time:", round(sum(batch_results$processing_time), 1), "seconds"),
    "",
    "DETAILED RESULTS",
    "----------------"
  )
  
  for (i in seq_len(nrow(batch_results))) {
    row <- batch_results[i, ]
    report_lines <- c(report_lines,
      paste0("\n[", i, "] ", row$model_file),
      paste("    Status:", row$status),
      paste("    Illuminated Pixels:", row$n_illuminated_pixels),
      paste("    Processing Time:", round(row$processing_time, 2), "s")
    )
    
    if (row$projection_file != "") {
      report_lines <- c(report_lines,
        paste("    Projection Output:", basename(row$projection_file)))
    }
    
    if (row$reflectance_file != "") {
      report_lines <- c(report_lines,
        paste("    Reflectance Output:", basename(row$reflectance_file)))
    }
  }
  
  # Save to file if requested
  if (!is.null(output_file)) {
    writeLines(report_lines, output_file)
    message("Report saved to: ", output_file)
  }
  
  # Print to console
  cat(paste(report_lines, collapse = "\n"), "\n")
  
  invisible(report_lines)
}
