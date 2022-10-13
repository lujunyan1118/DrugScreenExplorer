load("testData/expectedOutput/screenData.RData")

testScreen <- readScreen(rawDir =  "testData/rawData/",
                         plateAnnotationFile = "testData/plateAnno/plateAnno.tsv",
                         wellAnnotationFile = "testData/wellAnno/wellAnno.tsv",
                         rowRange = c(3,18), colRange = 2,
                         sep = "[;\t]",
                         negWell <- c("DMSO","PBS"), posWell = c())


testScreen_norm <- readScreen(rawDir =  "testData/rawData/",
                                      plateAnnotationFile = "testData/plateAnno/plateAnno.tsv",
                                      wellAnnotationFile = "testData/wellAnno/wellAnno.tsv",
                                      rowRange = c(3,18), colRange = 2,
                                      sep = "[;\t]",
                                      negWell <- c("DMSO","PBS"), posWell = c(),
                                      normalization = TRUE, method = "negatives")

testScreen_norm_discard <- readScreen(rawDir =  "testData/rawData/",
                         plateAnnotationFile = "testData/plateAnno/plateAnno.tsv",
                         wellAnnotationFile = "testData/wellAnno/wellAnno.tsv",
                         rowRange = c(3,18), colRange = 2,
                         sep = "[;\t]",
                         negWell <- c("DMSO","PBS"), posWell = c(),
                         normalization = TRUE, method = "negatives", discardLayer = 2)



test_that("Read the whole screen data", {
  expect_equal(testScreen, screenData)
})

test_that("Read the whole screen data, with normalization", {
  expect_equal(testScreen_norm, screenData_norm)
})


test_that("Read the whole screen data, with normalization, discarding 2 outer layers", {
  expect_equal(testScreen_norm_discard, screenData_norm_discard)
})
