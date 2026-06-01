library(Seurat)
library(ggplot2)
library(openxlsx)

# =========================
# Create directories
# =========================

create_dirs <- function(base_dir) {

  dir.create(base_dir, showWarnings = FALSE)

  dir.create(
    paste0(base_dir, "/plots"),
    showWarnings = FALSE
  )

  dir.create(
    paste0(base_dir, "/metadata"),
    showWarnings = FALSE
  )

  dir.create(
    paste0(base_dir, "/markers"),
    showWarnings = FALSE
  )

  dir.create(
    paste0(base_dir, "/objects"),
    showWarnings = FALSE
  )
}

# =========================
# Save JPEG plots
# =========================

save_plot_jpeg <- function(
  plot_obj,
  filename,
  width = 1800,
  height = 1400,
  res = 200
) {

  jpeg(
    filename,
    width = width,
    height = height,
    res = res
  )

  print(plot_obj)

  dev.off()
}
