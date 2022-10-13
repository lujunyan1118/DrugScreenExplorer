load("testData/expectedOutput/expected_values.RData")

test_that("Read one 384-well plate", {
  expect_identical(readPlate("testData/example_plates/plate1.txt",
                         rowRange = c(3,18), colRange =2, sep = ";", commaAsDecimal = FALSE ),
               plate1)
})


test_that("Read one 1536-well plate", {
  expect_identical(readPlate("testData/example_plates/plate2.txt",
                             rowRange = c(1,32), colRange =1, sep = "\t", commaAsDecimal = FALSE ),
                   plate2)
})

