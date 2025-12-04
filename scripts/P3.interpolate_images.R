args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <file> <path>")
}

file = args[1]
path = args[2]

library(raster)
library(sp)

source("./scripts/custom_functions.R")

# image
img <- raster::stack(file.path(path, file))

ext <- extent(img)
ras_dim <- dim(img)
rasTemp <- raster(nrows = ras_dim[1], ncol = ras_dim[2], ext = ext)

# mask
mask <- raster::raster("images/damsel_mask.tif")
crs(mask) = NA
NA_th <- sum(is.na(mask[]))

# interpolate if NA count higher than NA threshold
iteration = 0
downsampled = FALSE
NA_count <- sum(is.na(img[[1]][]))
while(NA_count > NA_th) {
  img <- sp::disaggregate(img, fact = c(15,10), method = "bilinear")
  img <- raster::resample(img, rasTemp, method = "ngb")
  img <- raster::mask(img, mask)
  img[img < 0] = 0  
  img[img > 255] = 255
  
  new_NA_count = sum(is.na(img[[1]][]))
  if (new_NA_count == NA_count) {
    if (new_NA_count > NA_th) {
      downsampled = TRUE
      img2 <- raster::resample(img, raster(nrows = 50, ncol = ras_dim[2] * 50 / ras_dim[1] , ext = ext), method = "ngb")
      img2 <- raster::resample(img2, rasTemp, method = "ngb")
      img[[1]][is.na(img[[1]])] = img2[[1]][is.na(img[[1]])]
      img[[2]][is.na(img[[2]])] = img2[[2]][is.na(img[[2]])]
      img[[3]][is.na(img[[3]])] = img2[[3]][is.na(img[[3]])]
      img <- raster::mask(img, mask)
      img[img < 0] = 0  
      img[img > 255] = 255
    }
    break
  }
  iteration = iteration + 1
  NA_count = new_NA_count
}

cat("\nNumber of iterations:", iteration)
cat("\nDownsampled:", downsampled)
cat("\n")

img <- raster::mask(img, mask)

output_path <- file.path("images/Transformed_images_interpolated", paste0(i, ".tif"))

raster::writeRaster(img, filename = output_path, overwrite = TRUE)
cat("\nProcessed image saved as:", output_path)

