# animalreflect

An R package for integrating 3D animal models with hyperspectral imaging data to extract projection surfaces and retrieve spectral reflectance information based on solar geometry.

## Overview

**animalreflect** enables researchers to:
- Load 3D animal models in OBJ format
- Load hyperspectral images in SPE format (300-2100 nm wavebands)
- Calculate solar position based on date, time, and location
- Extract illuminated projection surfaces of 3D models
- Retrieve hyperspectral reflectance data from projected surfaces
- Export projection surfaces as GeoTIFF files
- Compute mean reflectance values at specific wavebands

## Features

- **3D Model Support**: Load and manipulate OBJ format animal models
- **Hyperspectral Integration**: Process SPE format hyperspectral images with multiple wavebands
- **Solar Position Calculation**: Time-varying solar geometry calculations supporting date/time/location inputs
- **Projection Extraction**: Extract illuminated projection surfaces based on solar angles
- **Reflectance Analysis**: Extract spectral reflectance data and compute statistics
- **Spatial Output**: Export results as GeoTIFF format with proper spatial referencing
- **Batch Processing**: Process multiple animal models efficiently

## Installation

```r
# Install from GitHub
devtools::install_github("mayinsitan1/animalreflect")
```

## Quick Start

```r
library(animalreflect)

# 1. Load 3D model and hyperspectral image
model <- load_obj_model("path/to/animal.obj")
hyper_data <- load_hyperspectral_spe("path/to/image.spe")

# 2. Calculate solar position (e.g., for a specific date and time)
solar <- calculate_solar_position(
  date = "2026-06-01",
  time = "12:00:00",
  latitude = 40.7128,
  longitude = -74.0060,
  azimuth = 180  # South-facing
)

# 3. Extract projection surface
projection <- extract_projection_surface(
  model = model,
  solar_position = solar
)

# 4. Extract reflectance data from projection surface
reflectance_data <- extract_reflectance_surface(
  projection = projection,
  hyperspectral_data = hyper_data,
  wavelengths = hyper_data$wavelengths
)

# 5. Export projection as GeoTIFF
export_projection_tif(
  projection = projection,
  output_file = "projection_output.tif"
)

# 6. Get mean reflectance at specific wavebands
mean_refl <- get_mean_reflectance(
  reflectance_data = reflectance_data,
  wavelengths = c(550, 850, 1650),  # RGB, NIR, SWIR
  hyperspectral_wavelengths = hyper_data$wavelengths
)
```

## Data Format Requirements

### OBJ Model Format
- Standard Wavefront OBJ format
- Can include vertex normals for improved projection calculations
- Should contain closed surface geometry

### SPE Format (Hyperspectral Images)
- Binary format with hyperspectral data
- Wavebands: typically 300-2100 nm range
- Spatial resolution: < 1 mm
- Metadata should include wavelength information

### Solar Angle Input
- **Azimuth**: Direction (0-360°, where 0° = North, 90° = East, 180° = South, 270° = West)
- **Elevation**: Sun angle above horizon (0-90°)

## Key Functions

- `load_obj_model()`: Load 3D OBJ model files
- `load_hyperspectral_spe()`: Load SPE format hyperspectral images
- `calculate_solar_position()`: Calculate solar position from date/time/location
- `extract_projection_surface()`: Extract illuminated projection surface
- `extract_reflectance_surface()`: Map hyperspectral data to projection surface
- `get_mean_reflectance()`: Calculate mean reflectance at specific wavelengths
- `export_projection_tif()`: Export projection as GeoTIFF
- `batch_process_animals()`: Process multiple models in batch mode

## Output Formats

- **Projection Surface (TIF)**: GeoTIFF format with spatial referencing
- **Reflectance Data**: RData format preserving full spectral information
- **Statistics**: CSV files with mean reflectance values

## System Requirements

- R ≥ 4.1
- Compiled dependencies for geometric computations

## Package Dependencies

Core dependencies include:
- `rgl`: 3D visualization and model manipulation
- `terra`/`raster`: Raster data handling
- `geometry`: Computational geometry operations
- `Rcpp`/`RcppArmadillo`: High-performance C++ computations
- `lubridate`: Date/time handling
- `sf`: Spatial data framework

## Contributing

Contributions are welcome! Please submit issues and pull requests to the [GitHub repository](https://github.com/mayinsitan1/animalreflect).

## License

MIT License - see LICENSE file for details

## Author

mayinsitan1

## Citation

If you use this package in your research, please cite it as:

```bibtex
@software{animalreflect2026,
  title = {animalreflect: Integration of 3D Animal Models with Hyperspectral Imaging},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/mayinsitan1/animalreflect}
}
```

## References

- OBJ File Format: https://en.wikipedia.org/wiki/Wavefront_.obj_file
- Hyperspectral Imaging: Hagen, N., & Kudenov, M. W. (2013). Review of snapshot spectral imaging technologies. Optical Engineering, 52(9), 090901.
- Solar Position Algorithm: Reda, I., & Andreas, A. (2008). Solar position algorithm for solar radiation applications. National Renewable Energy Laboratory, 1-34.