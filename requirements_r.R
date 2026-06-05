packages <- c("dplyr", "data.table", "twingp", "readr", "tidyr", "hetGP")

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}