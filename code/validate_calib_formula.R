###this is code is for validating
## beta r = 2a Es^2/sqrt(N)

N <- 1e5
s <- rnorm(N,mean = 0, sd = 2)
a <- 0.07
mu <- a*s+1/2
print(max(abs(mu)))
mu[mu > 1] <- 1
mu[mu < 0] <- 0
epsilon <- rbinom(N,1,mu)
epsilon <- as.numeric(epsilon == 1)*(1-mu) + as.numeric(epsilon == 0)*(-mu)
y <- a*s + epsilon
y <- y/(sd(y)*sqrt(N))
mean(y)
var(y)
mean(y*s)
2*a*2^2/sqrt(N)
