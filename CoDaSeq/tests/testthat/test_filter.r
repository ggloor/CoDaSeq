library(CoDaSeq)
data(ak_op)

context("core functions")

filt <- codaSeq.filter(ak_op, min.reads=1000, min.prop=0.01, max.prop=1,
    min.occurrence=0.25, samples.by.row=FALSE)

f.n0 <- cmultRepl(t(filt), label=0, method="CZM")

f.clr <- codaSeq.clr(f.n0, samples.by.row=TRUE)


test_that("filter function is sane", {
  expect_equal(nrow(filt), 167)
  expect_equal(ncol(filt), 30)
  expect_equal(min(filt), 0)

})

test_that("clr function is sane", {
  expect_equal(ncol(f.clr), 167)
  expect_equal(nrow(f.clr), 30)
  expect_equal(anyNA(f.clr), FALSE)

})

