context("test input checking for mipred.cv")

cll_bin <- readRDS("CLL_BMJ.rds")

cll_bin$srv5y_s[cll_bin$srv5y>12] <- 0  # Apply an administrative censorship at t=84 months
cll_bin$srv5y[cll_bin$srv5y>12]  <- 12

cll_bin$Status[cll_bin$srv5y_s==1]<- 1  # Define the "Status" variable
cll_bin$Status[cll_bin$srv5y_s==0] <- 0  # Ss numeric -> 1:Dead, 0:Alive

cll_bin$Censor <- NULL
cll_bin$srv5y <- NULL
cll_bin$srv5y_s <- NULL

test_that("formula not recognized or input",{
    expect_error(mipred.cv(data=cll_bin[1:500,-1], nimp=5, folds=5), "formula not provided")
    expect_error(mipred.cv(12345,data=cll_bin[1:500,-1], nimp=5, folds=5), "formula argument not recognized as formula")
})

family<-NULL
test_that("family not recognized or input",{
  expect_error(mipred.cv(Status~age10, data=cll_bin[1:500,-1], nimp=5, folds=5), "argument \"family\" is missing, with no default")
  expect_error(mipred.cv(Status~age10, family, data=cll_bin[1:500,-1], nimp=5, folds=5), "family not recognized")
})

data1 <- matrix(5,nrow=2,ncol=2)
data2 <- cll_bin
names(data2)[3]<-"age10"

test_that("data properly input",{
  expect_error(mipred.cv(Status~age10, family=binomial,  nimp=5, folds=5), "data not provided")
  expect_error(mipred.cv(Status~age10, family=binomial, data=data1, nimp=5, folds=5), "data should be a data frame")
  expect_error(mipred.cv(Status~age10, family=binomial, data=data2, nimp=5, folds=5), "Duplicate names found in data: age10")
})

test_that("nimp properly input",{
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], folds=5), "nimp not provided")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=-5, folds=5), "nimp must be integer number greater than zero")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp="a", folds=5), "nimp must be integer number greater than zero")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=c(1,2,5), folds=5), "nimp must be integer number greater than zero")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1.24, folds=5), "nimp must be integer number greater than zero")
})

test_that("folds properly input",{
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds=c(1,5)), "folds must be integer number greater than 1")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds=-5), "folds must be integer number greater than 1")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds="a"), "folds must be integer number greater than 1")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds=0), "folds must be integer number greater than 1")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds=1.45), "folds must be integer number greater than 1")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=5, folds=3000), "folds too large")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=1, mice.options=mice.options), "folds must be integer number greater than 1")
})

mice.options <- NULL
mice.options$m <- 5
mice.options$maxit=1
test_that("argument m in mice.options is ignored",{
  expect_warning(mipred.cv(Status~age10+cyto, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=2, mice.options=mice.options), "mice.options argument m ignored, using nimp instead")
})

test_that("proper specification of seed in argument mice.options",{
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=2, mice.options=list(maxit=1, seed=matrix(c(1,2)))), "in mice.options seed must be a vector")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=2, mice.options=list(maxit=1, seed="a")), "in mice.options seed must be a numeric vector")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=2, mice.options=list(maxit=1, seed=c(1.22,3.44))), "in mice.options seed must be positive integer")
  expect_error(mipred.cv(Status~age10, family=binomial, data=cll_bin[1:500,-1], nimp=1, folds=2, mice.options=list(maxit=1, seed=-3)), "in mice.options seed must be positive integer")
})





