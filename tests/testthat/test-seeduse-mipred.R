context("use of seed")

skip('seeduse tests mipred')

cll <- readRDS("CLL_BMJ.rds")

cll_bin<-cll

cll_bin$srv5y_s[cll_bin$srv5y>12] <- 0  # Apply an administrative censorship at t=84 months
cll_bin$srv5y[cll_bin$srv5y>12]  <- 12

cll_bin$Status[cll_bin$srv5y_s==1]<- 1  # Define the "Status" variable
cll_bin$Status[cll_bin$srv5y_s==0] <- 0  # Ss numeric -> 1:Dead, 0:Alive

cll_bin$Censor <- NULL
cll_bin$srv5y <- NULL
cll_bin$srv5y_s <- NULL

seeds <- c(7350, 5880, 2571, 3887, 9745, 5608, 3468, 2711, 3448, 4676,
  1627, 5087, 6479, 3901, 1037, 6622, 4832, 3382, 8515, 268, 3210,
  5944, 3380, 3646, 9188, 2750, 7987, 466, 2981, 136, 5134, 8476,
  5618, 8948, 1590, 5281, 2551, 2052, 1444, 8688, 1691, 6398, 4728,
  8267, 5274, 4154, 982, 7540, 3081, 5286, 6405, 1238, 5589, 8221,
  6026, 1922, 7651, 3100, 4341, 2134, 9685, 4524, 9750, 7618, 8211,
  2438, 824, 6982, 6548, 9743, 4302, 9038, 185, 4418, 1306, 6643,
  9865, 3436, 6151, 9647, 4261, 8109, 7573, 127, 571, 4545, 7245,
  4764, 5038, 9037, 6085, 1459, 5124, 9463, 942, 2241, 6492, 3190,
  8229, 7099)

output <- list(pred = structure(c(0.272822624187677, 0.26335349393881,
  0.245055602778557, 0.316090037563377, 0.313007499662245, 0.297912651367598,
  0.307891918762763, 0.255163876316593, 0.331122267440826, 0.364964866709756
), .Dim = c(5L, 2L), .Dimnames = list(c("501", "502", "503",
  "504", "505"), NULL)), linpred = structure(c(-0.980348575240286,
    -1.02861092500769, -1.1251589655607, -0.771799173093083, -0.786096224788477,
    -0.857257497214206, -0.809993321899551, -1.07125823500496, -0.703113586720463,
    -0.553879232018371), .Dim = c(5L, 2L), .Dimnames = list(c("501",
      "502", "503", "504", "505"), NULL)))

mice.options <- NULL
mice.options$seed <- as.vector(matrix(seeds,ncol=1,nrow=5*5))[1:4]
mice.options$maxit <- 2

set.seed(12345)
test_that("output replicates from specified seed for mipred and averaging method",{
  expect_equal(mipred(Status~age10+cyto,family=binomial,data=cll_bin[1:500,-1],newdata=cll_bin[501:505,c(-1,-10)], nimp=2, folds=2, mice.options=mice.options)[2:3],output)
    })

output <- list(pred = structure(c(0.313619237134538, 0.287549302953276,
  0.243573611866073, 0.313856654363961, 0.278459046089519, 0.313619237134538,
  0.287549302953276, 0.243573611866073, 0.313856654363961, 0.324948659032894
), .Dim = c(5L, 2L), .Dimnames = list(c("501", "502", "503",
  "504", "505"), NULL)), linpred = structure(c(-0.783252892154253,
    -0.907316378553107, -1.13318601836123, -0.782150197640975, -0.952118141968731,
    -0.783252892154253, -0.907316378553107, -1.13318601836123, -0.782150197640975,
    -0.731121551027735), .Dim = c(5L, 2L), .Dimnames = list(c("501",
      "502", "503", "504", "505"), NULL)))

mice.options <- NULL
mice.options$seed <- as.vector(matrix(seeds,ncol=1,nrow=5*5))[1:2]
mice.options$maxit <- 2

set.seed(12345)
test_that("output replicates from specified seed for mipred and rubin method",{
  expect_equal(mipred(Status~age10+cyto,family=binomial,data=cll_bin[1:500,-1],newdata=cll_bin[501:505,c(-1,-10)], nimp=2, folds=2, mice.options=mice.options, method="rubin")[2:3],output)
})
