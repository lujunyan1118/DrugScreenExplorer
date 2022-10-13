load("testData/expectedOutput/screenData.RData")
load("testData/expectedOutput/edgeFactor.RData")

testEdgeFactor_loess <- fitEdgeEffect(screenData_norm, method = "loess",useNeg = TRUE, useLowConcentrations = 0)$edgeFactor
testEdgeFactor_loess_allWell <- fitEdgeEffect(screenData_norm, method = "loess",useNeg = FALSE, useLowConcentrations = 0)$edgeFactor
testEdgeFactor_sigmoid <- fitEdgeEffect(screenData_norm, method = "sigmoid",useNeg = TRUE, useLowConcentrations = 0)$edgeFactor
testEdgeFactor_sigmoid_allWell <- fitEdgeEffect(screenData_norm, method = "sigmoid",useNeg = FALSE, useLowConcentrations = 0)$edgeFactor


test_that("Fit edge factor using loess method with negative control wells", {
  expect_equal(edgeFactor_loess, testEdgeFactor_loess)
})

test_that("Fit edge factor using loess method with all wells", {
  expect_equal(edgeFactor_loess_allWell, testEdgeFactor_loess_allWell)
})

test_that("Fit edge factor using sigmoid method with negative control wells", {
  expect_equal(edgeFactor_sigmoid, testEdgeFactor_sigmoid)
})

test_that("Fit edge factor using loess method with all wells", {
  expect_equal(edgeFactor_sigmoid_allWell, testEdgeFactor_sigmoid_allWell)
})


test_that("Estimate edge effect strength using 1 layer", {
  expect_equal(estimateOnePlate(screenOnePlate, nLayer =1, onlyNeg = FALSE), edgeVal1)
})

test_that("Estimate edge effect strength using 2 layers and only negative wells", {
  expect_equal(estimateOnePlate(screenOnePlate, nLayer =2, onlyNeg = TRUE), edgeVal2)
})
