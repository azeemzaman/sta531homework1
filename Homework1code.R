library(BB)
library(truncnorm)
f <- function(y, theta){
  return((y-theta)/(1+(y-theta)^2))
}

g <- function(theta){
  y.vec <- c(-2,-1,0,1.5,2.5)
  return(sum(sapply(y.vec, function(y) {f(y,theta)})))
}

pos.mod <- BBoptim(1, g)$value

d2.log.f <- function(y, theta){
  return(((y-theta)^2-1)/((y-theta)^2+1)^2)
}

d2.log.g <- function(y.vec, theta){
  return(sum(sapply(y.vec, function(y) {d2.log.f(y,theta)})))
}
y.vec <- c(-2,-1,0,1.5,2.5)
nor.var <- -1/d2.log.g(y.vec, pos.mod)
# dist is N(pos.mod, nor.var)

# plot actual density
seq.length = 10000
theta.seq <- seq(from = 0,
                 to = 1,
                 length.out = seq.length)
cauchy.f <- function(y, theta){
  return((1+(y-theta)^2)^(-2))
}
unnorm.post <- function(y.vec, theta){
  return(prod(sapply(y.vec, function(y) {cauchy.f(y, theta)})))
}

unnorm.den <- sapply(theta.seq, function(theta) {unnorm.post(y.vec, theta)})
plot(theta.seq, unnorm.den)
post.tot <- sum(unnorm.den/seq.length)
norm.const <- 1/post.tot
norm.post <- function(y.vec, theta){
  norm.const*unnorm.post(y.vec, theta)
}
norm.den <- sapply(theta.seq, function(theta) {norm.post(y.vec, theta)})
plot(theta.seq, norm.den, type = "l")
lines(theta.seq, dtruncnorm(theta.seq, 
                      a = 0,
                      b = 1,
                      mean = pos.mod, 
                      sd = sqrt(nor.var)))
