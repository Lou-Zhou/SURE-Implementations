library(testthat)
library(tidyverse)
source("implementations.R")
check_equal_pca <- function(X, tolerance = 6) {
  check_rotation = all(round(run_pca(X)$rotation, tolerance) == 
        round(prcomp(X, center = TRUE, scale. = TRUE)$rotation, tolerance))
  check_sdev = all(round(run_pca(X)$sdev, tolerance) == 
                     round(prcomp(X, center = TRUE, scale. = TRUE)$sdev, tolerance))
  check_x = all(round(run_pca(X)$x, tolerance) == 
                  round(prcomp(X, center = TRUE, scale. = TRUE)$x, tolerance))
  check_rotation & check_sdev & check_x
}
test_that("PCA random matrix", {
  r <- sample(1:1000, 2, replace = TRUE) # Generate random size
  random_matrix <- matrix(runif(r[1] * r[2]), nrow = r[1])
  expect_true(check_equal_pca(random_matrix))
})

test_that("Starbucks PCA", {
  starbucks <- read_csv(
    "https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2021/2021-12-21/starbucks.csv"
  ) |>
    mutate(trans_fat_g = as.numeric(trans_fat_g), fiber_g = as.numeric(fiber_g))  |>
    select(serv_size_m_l:caffeine_mg)
   expect_true(check_equal_pca(data.matrix(starbucks)))
})
