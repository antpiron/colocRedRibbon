context("RedRibbonColoc")

library(data.table)
library(magrittr)

n <- 100
a <- runif(n)
b <- runif(n)
lead <- min(a, b) * 1e-16
frac <- floor(0.1*n)
ld <-   (frac:1) / frac
a[1:frac] <- exp(ld * log(lead) + (1-ld) * log(a[1:frac]))
b[1:frac] <- exp(ld * log(lead) + (1-ld) * log(b[1:frac]))

data <- data.table(id=paste0("g", 1:n),
                   a=a, b=b,
                   position=(1:n) * 100)

rrColoc <- RedRibbonColoc(data)

test_that("RedRibbon", {
    expect_equal(rrColoc$quadrants$whole$positions[1], 1)
    expect(rrColoc$quadrants$whole$log_padj >= -log(0.05), "padj > 0.05")
})

data[, a.eaf := c(0.5, runif(n-1))]
data[, b.eaf := a.eaf]
data[, a.n := 404]
data[, b.n := 404]

rrColoc <- RedRibbonColoc(data) %>% coloc()

test_that("Coloc", {
    expect_equal(rrColoc$coloc$bestSnp, "g1")
    expect(rrColoc$coloc$PP.H4.abf >= 0.8, "PP.H4.abf < 0.8")
    expect(rrColoc$coloc$SNP.PP.H4 >= 0.8, "SNP.PP.H4 < 0.8")
})

data[, a.dir := c(1,  1, rnorm(n-2))]
data[, b.dir := c(1, -1, rnorm(n-2))]

rrColoc <- RedRibbonColoc(data, risk="a") %>% coloc()
test_that("Coloc", {
    expect_equal(rrColoc$coloc$bestSnp, "g1")
    expect(rrColoc$coloc$PP.H4.abf >= 0.8, "PP.H4.abf < 0.8")
    expect(rrColoc$coloc$SNP.PP.H4 >= 0.8, "SNP.PP.H4 < 0.8")
})

rrColoc <- RedRibbonColoc(data, risk="a", effect=`<=`) %>% coloc()
test_that("Coloc", {
    expect_equal(rrColoc$coloc$bestSnp, "g2")
    expect(rrColoc$coloc$PP.H4.abf >= 0.8, "PP.H4.abf < 0.8")
    expect(rrColoc$coloc$SNP.PP.H4 >= 0.8, "SNP.PP.H4 < 0.8")
})
