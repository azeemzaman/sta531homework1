# Problem 9 code#

# MLE function
MLE <- function(y){
  if ((y <= 1) && (y >= 0)){
    return(y)
  } else if (y < 0){
    return(0)
  } else {
    return(1)
  }
}

postMean <- function(y){
  post.mean <- y + (dnorm(y)-dnorm(1-y))/(pnorm(1-y)-pnorm(-y))
  return(post.mean)
}

true.theta <- runif(1)
trails <- 10000
variance <- seq(from = .1, to = 5, length.out = 100)
computeMSE <- function(true.theta, 
                       trails, 
                       variance,
                       to.return = "MSE"){
  thetas <- rnorm(trails,
                  mean = true.theta,
                  sd = sqrt(variance))
  MLE.res <- sapply(thetas, MLE)
  post.mean.res <- sapply(thetas, postMean)
  MLE.mean <- mean(MLE.res)
  post.mean.mean <- mean(post.mean.res)
  MLE.var <- var(MLE.res)
  post.mean.var <- var(post.mean.res)
  MLE.MSE <- (MLE.mean - true.theta)^2+MLE.var
  post.mean.MSE <- ((post.mean.mean - true.theta)^2+
                    post.mean.var)
  MLE.vals <- list("results" = MLE.res,
                   "mean" = MLE.mean,
                   "variance" = MLE.var,
                   "MSE" = MLE.MSE)
  post.mean.vals <- list("results" = post.mean.res, 
                         "mean" = post.mean.mean,
                     "variance" = post.mean.var, 
                     "MSE" = post.mean.MSE)
  return(c(MLE.vals[[to.return]], 
           post.mean.vals[[to.return]]))
}
makeLists <- function(true.theta,
                      trails,
                      variance, 
                      to.return){
  list.of.vals <- lapply(variance, function(x){
    computeMSE(true.theta, trails, x, to.return)
  })
  MLE.vals <- sapply(list.of.vals, "[[", 1)
  post.mean.vals <- sapply(list.of.vals, "[[", 2)
  matplot(variance, 
       cbind(MLE.vals, post.mean.vals),
       type = "l",
       lty = c(1,2),
       main = paste0("Plot of ", to.return),
       ylab = to.return)
}

makeLists(true.theta, trails, variance, "MSE")
makeLists(1, trails, variance, "MSE")
makeLists(0, trails, variance, "MSE")
makeLists(true.theta, trails, variance, "variance")

# Problem 15
c <- qnorm(.75)
val1 <- (5/4)*(1+sqrt(5)*c/2)
val2 <- (5/4)*(1-sqrt(5)*c/2)
p1 <- 1-pnorm(val1, mean = 1, sd = 1)
p2 <- pnorm(val2, mean = 1, sd = 1)
prob <- p1+p2

computeCoverage <- function(cover, true.theta){
  c <- qnorm((1+cover)/2)
  val1 <- (5/4)*(true.theta + sqrt(5)*c/2)
  val2 <- (5/4)*(true.theta - sqrt(5)*c/2)
  p1 <- 1-pnorm(val1, mean = true.theta, sd = 1)
  p2 <- pnorm(val2, mean = true.theta, sd = 1)
  prob <- p1+p2
  return(prob)
}
theta.seq <- seq(from = -16, to = 16, length.out = 10000)
coverage <- sapply(theta.seq, function(x) {
  computeCoverage(.5, x)})
plot(theta.seq, coverage, type = "l")
hist(coverage)
computeCoverage(.5, 2)

computeFreq <- function(cover, true.theta){
  c <- qnorm((1+cover)/2)
  val1 <- true.theta+c
  val2 <- true.theta-c
  p1 <- 1- pnorm(val1, mean = true.theta, sd = 1)
  p2 <- pnoom(val2, mean = true.theta, sd = 1)
  prob <- p1+p2
  return(prob)
}
freq.cover <- sapply(theta.seq, function(x) {
  computeFreq(.5, x)
})
lines(theta.seq, freq.cover)