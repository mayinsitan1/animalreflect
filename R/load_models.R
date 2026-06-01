#' Load 3D Animal Model from OBJ File
#'
#' Loads a Wavefront OBJ file containing a 3D animal model.
#'
#' @param obj_file Character string specifying the path to the OBJ file.
#' @param scale Numeric. Optional scaling factor to apply to the model coordinates.
#'   Default is 1 (no scaling).
#'
#' @details
#' This function reads an OBJ file and extracts vertices, faces, and normals.
#' The resulting object can be used with other functions in the package for
#' projection extraction and analysis.
#'
#' @return A list with class "obj_model" containing:
#'   \item{vertices}{Numeric matrix of vertex coordinates (x, y, z)}
#'   \item{faces}{Integer matrix of face definitions (triangles)}
#'   \item{normals}{Numeric matrix of vertex normals (if available)}
#'   \item{filename}{Original filename}
#'   \item{n_vertices}{Number of vertices}
#'   \item{n_faces}{Number of faces}
#'   \item{scale_factor}{Applied scaling factor}
#'
#' @examples
#' \dontrun{
#'   model <- load_obj_model("path/to/animal.obj")
#'   print(model)
#'   summary(model)
#' }
#'
#' @export
#' @importFrom utils read.table
load_obj_model <- function(obj_file, scale = 1) {
  
  if (!file.exists(obj_file)) {
    stop("OBJ file not found: ", obj_file)
  }
  
  # Read OBJ file
  lines <- readLines(obj_file)
  
  vertices <- NULL
  faces <- NULL
  normals <- NULL
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Skip comments and empty lines
    if (grepl("^#", line) || line == "") next
    
    parts <- strsplit(line, "\\s+")[[1]]
    
    # Parse vertices
    if (parts[1] == "v" && length(parts) >= 4) {
      v <- as.numeric(parts[2:4]) * scale
      vertices <- rbind(vertices, v)
    }
    
    # Parse vertex normals
    else if (parts[1] == "vn" && length(parts) >= 4) {
      vn <- as.numeric(parts[2:4])
      normals <- rbind(normals, vn)
    }
    
    # Parse faces
    else if (parts[1] == "f" && length(parts) >= 4) {
      # Handle different face formats (v, v/vt, v/vt/vn, v//vn)
      face_verts <- integer()
      
      for (j in 2:length(parts)) {
        indices <- strsplit(parts[j], "/")[[1]]
        face_verts <- c(face_verts, as.integer(indices[1]))
      }
      
      # Convert to triangles if face has more than 3 vertices
      if (length(face_verts) == 3) {
        faces <- rbind(faces, face_verts)
      } else if (length(face_verts) > 3) {
        # Fan triangulation
        for (k in 2:(length(face_verts) - 1)) {
          faces <- rbind(faces, c(face_verts[1], face_verts[k], face_verts[k + 1]))
        }
      }
    }
  }
  
  # Convert to numeric matrices if NULL
  if (is.null(vertices)) {
    stop("No vertices found in OBJ file")
  }
  if (is.null(faces)) {
    stop("No faces found in OBJ file")
  }
  
  # Create object
  obj_model <- list(
    vertices = vertices,
    faces = faces,
    normals = normals,
    filename = basename(obj_file),
    n_vertices = nrow(vertices),
    n_faces = nrow(faces),
    scale_factor = scale
  )
  
  class(obj_model) <- "obj_model"
  
  return(obj_model)
}


#' Load Hyperspectral Image from SPE File
#'
#' Loads a hyperspectral image stored in SPE (Princeton Instruments) format.
#'
#' @param spe_file Character string specifying the path to the SPE file.
#' @param wavelength_file Optional character string specifying path to a file
#'   containing wavelength information. If NULL, wavelengths are extracted from
#'   SPE metadata if available.
#'
#' @details
#' SPE files are binary format files used by Princeton Instruments spectrometers.
#' This function reads the binary data and extracts the hyperspectral cube and
#' associated metadata.
#'
#' The expected data structure is:
#' - Spatial dimensions: X (horizontal), Y (vertical)
#' - Spectral dimension: wavelength/wavenumber bands (300-2100 nm typical)
#'
#' @return A list with class "hyperspectral_spe" containing:
#'   \item{data}{3D numeric array with dimensions (x, y, bands)}
#'   \item{wavelengths}{Numeric vector of wavelength values in nanometers}
#'   \item{x_pixels}{Number of pixels in X dimension}
#'   \item{y_pixels}{Number of pixels in Y dimension}
#'   \item{n_bands}{Number of spectral bands}
#'   \item{metadata}{List of metadata from SPE file}
#'   \item{filename}{Original filename}
#'
#' @examples
#' \dontrun{
#'   hyper_data <- load_hyperspectral_spe("path/to/image.spe")
#'   str(hyper_data)
#'   summary(hyper_data$wavelengths)
#' }
#'
#' @export
load_hyperspectral_spe <- function(spe_file, wavelength_file = NULL) {
  
  if (!file.exists(spe_file)) {
    stop("SPE file not found: ", spe_file)
  }
  
  # Read binary SPE file
  con <- file(spe_file, "rb")
  on.exit(close(con))
  
  # SPE file structure (simplified version for common layouts)
  # Typically: header (4100 bytes) + footer + data
  
  # Read header to get dimensions
  header <- readBin(con, what = "raw", n = 4100)
  
  # Extract relevant header information
  # Note: This is a simplified parser. SPE format is complex.
  # Real implementation should use external libraries or detailed format specs
  
  # Read full file content
  seek(con, 0)
  full_data <- readBin(con, what = "raw", n = file.size(spe_file))
  
  # Placeholder: assume data structure needs proper SPE format parser
  # For now, return structure with message
  message("Note: SPE file parsing requires proper SPE format documentation.")
  message("Current implementation is a placeholder.")
  
  # Create a basic structure (to be filled with actual SPE parsing)
  hyper_object <- list(
    data = NULL,  # 3D array to be filled
    wavelengths = NULL,
    x_pixels = NA,
    y_pixels = NA,
    n_bands = NA,
    metadata = list(filename = spe_file),
    filename = basename(spe_file)
  )
  
  class(hyper_object) <- "hyperspectral_spe"
  
  return(hyper_object)
}


#' Print method for obj_model
#'
#' @param x An obj_model object
#' @param ... Additional arguments passed to print
#'
#' @export
print.obj_model <- function(x, ...) {
  cat("3D Animal Model (OBJ Format)\n")
  cat("============================", "\n")
  cat("File:", x$filename, "\n")
  cat("Vertices:", x$n_vertices, "\n")
  cat("Faces:", x$n_faces, "\n")
  cat("Scale factor:", x$scale_factor, "\n")
  cat("Vertex coordinate range:\n")
  print(apply(x$vertices, 2, range))
}


#' Print method for hyperspectral_spe
#'
#' @param x A hyperspectral_spe object
#' @param ... Additional arguments passed to print
#'
#' @export
print.hyperspectral_spe <- function(x, ...) {
  cat("Hyperspectral Image (SPE Format)\n")
  cat("================================\n")
  cat("File:", x$filename, "\n")
  cat("Dimensions (X x Y x Bands):", 
      x$x_pixels, "x", x$y_pixels, "x", x$n_bands, "\n")
  if (!is.null(x$wavelengths)) {
    cat("Wavelength range:", min(x$wavelengths), "-", max(x$wavelengths), "nm\n")
  }
}
