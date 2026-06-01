#' Extract Illuminated Projection Surface
#'
#' Extracts the projection surface of a 3D animal model based on solar geometry.
#'
#' @param model An obj_model object loaded via \code{\link{load_obj_model}}.
#' @param solar_position A solar_position object from \code{\link{calculate_solar_position}}.
#' @param include_shadow Logical. If TRUE, also computes shadow projection.
#'   Default is FALSE.
#'
#' @details
#' This function projects the 3D model onto a 2D plane perpendicular to the
#' solar radiation direction. The projection is computed by:
#'
#' 1. Calculating the projection direction (inverse of solar vector)
#' 2. Finding all faces that are front-facing to the sun
#'   (normal dot product with solar vector > 0)
#' 3. Projecting vertices onto the plane perpendicular to solar direction
#' 4. Creating a 2D raster of the projected surface
#'
#' @return A list with class "projection_surface" containing:
#'   \item{projected_vertices}{Matrix of 2D projected vertex coordinates}
#'   \item{projected_faces}{Face indices for the projected surface}
#'   \item{illuminated_faces}{Logical vector indicating front-facing faces}
#'   \item{raster}{2D raster representation of the projection}
#'   \item{solar_position}{The solar_position object used}
#'   \item{metadata}{List with projection information}
#'
#' @examples
#' \dontrun{
#'   model <- load_obj_model("animal.obj")
#'   solar <- calculate_solar_position(
#'     date = "2026-06-01",
#'     time = "12:00:00",
#'     latitude = 40.7128,
#'     longitude = -74.0060
#'   )
#'   projection <- extract_projection_surface(model, solar)
#'   plot(projection)
#' }
#'
#' @export
extract_projection_surface <- function(model, solar_position, include_shadow = FALSE) {
  
  if (!inherits(model, "obj_model")) {
    stop("model must be an obj_model object")
  }
  if (!inherits(solar_position, "solar_position")) {
    stop("solar_position must be a solar_position object")
  }
  
  # Extract solar direction (inverse of solar vector)
  solar_vector <- solar_position$solar_vector
  projection_direction <- -solar_vector / sqrt(sum(solar_vector ^ 2))
  
  # Calculate face normals using cross product
  vertices <- model$vertices
  faces <- model$faces
  n_faces <- nrow(faces)
  face_normals <- matrix(0, nrow = n_faces, ncol = 3)
  
  for (i in seq_len(n_faces)) {
    v1 <- vertices[faces[i, 1], ]
    v2 <- vertices[faces[i, 2], ]
    v3 <- vertices[faces[i, 3], ]
    
    edge1 <- v2 - v1
    edge2 <- v3 - v1
    
    normal <- crossprod_3d(edge1, edge2)
    norm_length <- sqrt(sum(normal ^ 2))
    
    if (norm_length > 0) {
      normal <- normal / norm_length
    }
    
    face_normals[i, ] <- normal
  }
  
  # Determine front-facing faces (dot product > 0)
  illuminated <- rep(FALSE, n_faces)
  for (i in seq_len(n_faces)) {
    dot_product <- sum(face_normals[i, ] * solar_vector)
    illuminated[i] <- dot_product > 0
  }
  
  # Create orthonormal basis for projection plane
  # Use solar vector as primary direction
  primary <- projection_direction
  
  # Choose secondary vector perpendicular to primary
  if (abs(primary[1]) < 0.9) {
    secondary <- c(1, 0, 0)
  } else {
    secondary <- c(0, 1, 0)
  }
  
  # Gram-Schmidt orthogonalization
  secondary <- secondary - sum(secondary * primary) * primary
  secondary <- secondary / sqrt(sum(secondary ^ 2))
  
  # Tertiary vector
  tertiary <- crossprod_3d(primary, secondary)
  
  # Project vertices onto 2D plane
  center <- apply(vertices, 2, mean)
  centered_vertices <- sweep(vertices, 2, center)
  
  projected_2d <- matrix(0, nrow = nrow(vertices), ncol = 2)
  for (i in seq_len(nrow(vertices))) {
    v <- centered_vertices[i, ]
    projected_2d[i, 1] <- sum(v * secondary)
    projected_2d[i, 2] <- sum(v * tertiary)
  }
  
  # Create raster from projection
  x_range <- range(projected_2d[, 1])
  y_range <- range(projected_2d[, 2])
  
  # Pad ranges slightly
  x_pad <- (x_range[2] - x_range[1]) * 0.05
  y_pad <- (y_range[2] - y_range[1]) * 0.05
  x_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  # Default raster resolution (can be adjusted)
  raster_res <- 100
  x_seq <- seq(x_range[1], x_range[2], length.out = raster_res)
  y_seq <- seq(y_range[1], y_range[2], length.out = raster_res)
  
  # Create raster grid
  raster_grid <- matrix(0, nrow = raster_res, ncol = raster_res)
  
  # Mark illuminated pixels
  for (i in seq_len(n_faces)) {
    if (!illuminated[i]) next
    
    # Get projected triangle vertices
    v1 <- projected_2d[faces[i, 1], ]
    v2 <- projected_2d[faces[i, 2], ]
    v3 <- projected_2d[faces[i, 3], ]
    
    # Mark pixels inside triangle (simplified - check bounding box)
    x_min <- max(1, which.min(abs(x_seq - min(v1[1], v2[1], v3[1]))))
    x_max <- min(raster_res, which.min(abs(x_seq - max(v1[1], v2[1], v3[1]))))
    y_min <- max(1, which.min(abs(y_seq - min(v1[2], v2[2], v3[2]))))
    y_max <- min(raster_res, which.min(abs(y_seq - max(v1[2], v2[2], v3[2]))))
    
    raster_grid[y_min:y_max, x_min:x_max] <- 1
  }
  
  # Create result object
  projection <- list(
    projected_vertices = projected_2d,
    projected_faces = faces,
    illuminated_faces = illuminated,
    raster = raster_grid,
    solar_position = solar_position,
    metadata = list(
      n_illuminated_faces = sum(illuminated),
      n_total_faces = n_faces,
      raster_resolution = raster_res,
      x_range = x_range,
      y_range = y_range,
      projection_basis = list(
        primary = primary,
        secondary = secondary,
        tertiary = tertiary
      )
    )
  )
  
  class(projection) <- "projection_surface"
  
  return(projection)
}


#' Cross product for 3D vectors
#'
#' @param a Numeric vector of length 3
#' @param b Numeric vector of length 3
#'
#' @return Numeric vector of length 3 representing the cross product
#'
#' @keywords internal
crossprod_3d <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}


#' Print method for projection_surface
#'
#' @param x A projection_surface object
#' @param ... Additional arguments passed to print
#'
#' @export
print.projection_surface <- function(x, ...) {
  cat("Projection Surface\n")
  cat("==================\n")
  cat("Illuminated faces:", x$metadata$n_illuminated_faces, "/",
      x$metadata$n_total_faces, "\n")
  cat("Raster resolution:", x$metadata$raster_resolution, "x",
      x$metadata$raster_resolution, "\n")
  cat("X range:", round(x$metadata$x_range[1], 2), "to",
      round(x$metadata$x_range[2], 2), "\n")
  cat("Y range:", round(x$metadata$y_range[1], 2), "to",
      round(x$metadata$y_range[2], 2), "\n")
}


#' Plot method for projection_surface
#'
#' @param x A projection_surface object
#' @param y Unused
#' @param ... Additional arguments passed to image
#'
#' @export
plot.projection_surface <- function(x, y, ...) {
  image(x$raster, main = "Projected Illuminated Surface",
        xlab = "X (normalized)", ylab = "Y (normalized)",
        col = c("white", "black"), ...)
}
