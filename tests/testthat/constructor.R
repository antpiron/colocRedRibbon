context("RedRibbonColoc")

library(data.table)

n <- 100
a <- runif(n)
b <- runif(n)
lead <- min(a, b) * 1e-16
frac <- floor(0.1*n)
ld <-   (frac:1) / frac
a[1:frac] <- exp(ld * log(lead) + (1-ld) * log(a[1:frac]))
b[1:frac] <- exp(ld * log(lead) + (1-ld) * log(b[1:frac]))

data <- data.table(id=paste0("g", 1:n), a=a, b=b, position=(1:n) * 100)

rrColoc <- RedRibbonColoc(data)

test_that("lead", {
  expect_equal(rrColoc$quadrants$whole$positions[1], 1)
})
