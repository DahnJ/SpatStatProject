### Analysis of Diggle's 1995 paper ###

## 1: Intro
# Spatial clustering: nonuniform geopgraphical distribution
# Temporal clustering: more illness in winter
# But unless there's a contagious element or spatio-temporally localised increases, these two operate separately

# Comparison to Knox - counting pairs that are within [0,s]x[0,t] of each other
# Similar results, except here it's a function of (s,t)

## 2: Spatio-temporal processes
# K(s,t) = 2*pi*s^2*t for homoPoisson
# If space-time independence, then K(s,t)=K(s)K(t)
# SO-intensity function l_2(s,t) - density at distance (s,t)
# K(s,t) - cummulative number till distance (s,t)

## 3: Estimation of second-order properties
# Why do the estimates look this way? Why is there {n(n-1)}^-1 just once?
# Edge-effects less important in temporal scale - why?
# That's why they define such a crude correction
# w,v=1 => ~Knox

## 4: Diagnostic for space-time clustering
# D = K - K1*K2
# D0 = D/(K1*K2) = K/(K1*K2) - 1
# Plot D to get information on the nature of spatio-temporal dependencies
# D0 has a nice interpretation as a proportional risk

# Third diagnostic
# D/sqrt(V) against K1*K2, similar to standardized residuals against actual values
# Two dimensional, which is nice, but the spatio-temporal scales are no longer explicit


## 5: Tests for space-time interactions
# Use Q(s,t) which is basically hat K and also Knox with weights - correction for edge effects
# U = \sum R(s,t)
# MC test for U
# Assumption of stationarity and radom permutations allow for non-stationarity


## 6: Examples
setwd("C:/gibbs/Baltimore")
library(splancs)
library(spatstat)

# Function for automatic plots. See ?stmctest for explanation of parameters
DigglePlots <- function(pts, times, poly, tlimits = c(0,1), s, tt, nsim = 100, quiet = FALSE){
  if(missing(poly)){ poly <- bboxx(spoints(c(0,1,0,1)))}
  if(missing(s)){ s <- seq(0,0.7,0.01)[-1]}
  if(missing(tt)) { tt <- s }
  
  stkh <- stkhat(pts, times, poly, tlimits, s, tt)
  stse <- stsecal(pts, times, poly, tlimits, s, tt)
  stmc <- stmctest(pts, times, poly, tlimits, s, tt, nsim, quiet)
  
  stdiagn(pts, stkh, stse, stmc)
  return(stkh)
}


## Poisson
lambda <- 100
poisson <- rpoispp3(lambda)
plot(poisson, main = "Poisson", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)
distances <- seq(0,0.7,0.01)[-1]
times <- distances
poisson.pts <- as.points(list(x = poisson$data$x, y = poisson$data$y))
poisson.times <- poisson$data$z

DigglePlots(poisson.pts, poisson.times)
paste("Poisson",toString(lambda),toString(nrow(poisson.pts)),sep="_")
paste("Realization of a Poisson point process. Intensity ", str(lambda), ". Number of points: ", toString(nrow(poisson.pts)), sep="")
paste("Diagnostic plots for a realization of a Poisson point process. Intensity ", str(lambda), ". Number of points: ", toString(nrow(poisson.pts)), sep="")



## Multidim Thomas process
lambdap <- 40
pois <- rpoispp3(lambdap)
points <- matrix(c(pois$data$x,pois$data$y,pois$data$z), ncol=3)

# Cropping to [0,1]^3
insideInt <- function(x,lower=0,upper=1){
  all(sapply(x,(function(y) (y >= lower) & (y <= upper))))
}

# Generating children for both 3-dim and "causal" Thomas processes
generateChildren <- function(p,lambda = 10, sigma = 0.01){
  n <- rpois(1,lambda)
  pts.uncropped <- t(t(matrix(rnorm(3*n,0,sigma),ncol=3)) + p)
  pts <- pts.uncropped[apply(pts.uncropped,1,insideInt),]
  return(pts)
}


lambdac <- 20
sigmac <- 0.05

Thomas <- apply(points,1,generateChildren, lambda = lambdac, sigma = sigmac)
Thomas <- do.call(rbind,Thomas)

Thomas.pts <- as.points(list(x = Thomas[,1], y = Thomas[,2]))
Thomas.times <- Thomas[,3]

plot(pp3(Thomas.pts[,1],Thomas.pts[,2],Thomas.times,as.box3(c(0,1,0,1,0,1))), 
     main = "Thomas", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

DigglePlots(Thomas.pts, Thomas.times)

paste("Thomas",toString(lambdap),toString(lambdac),gsub("\\.","p",toString(sigmac)),toString(nrow(Thomas.pts)),sep="_")
paste("Realization of a Thomas point process. Intensity of parent process ", toString(lambdap), ". Child process: intensity ", lambdap,", distribution \\$N(0,", toString(sigmac),")\\$."," Number of points: ", toString(nrow(Thomas.pts)), ".", sep="")
paste("Diagnostic plots for a realization of a Thomas point process. Intensity of parent process ", toString(lambdap), ". Child process: intensity ", lambdap,", distribution \\$N(0,", toString(sigmac),")\\$."," Number of points: ", toString(nrow(Thomas.pts)), ".", sep="")





## 3-dim Matern hard-core II
lambda <- 1000
pois <- rpoispp3(lambda)
points <- matrix(c(pois$data$x,pois$data$y,pois$data$z), ncol=3)
n <- nrow(points)

distances <- as.matrix(dist(points))
diag(distances) <- rep(100,n)

r <- 0.05
keep <- rep(TRUE,n)
for (i in 1:n){
  if (min(distances[i,]) < r) {
    # print(min(distances[i,]))
    keep[i] <- FALSE
    distances[i,] <- 100
    distances[,i] <- 100
  }
}
points <- points[keep,]


plot(pp3(points[,1],points[,2],points[,3],as.box3(c(0,1,0,1,0,1))), 
     main = "Matern", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

Matern.pts <- as.points(list(x=points[,1],y=points[,2]))
Matern.times <- points[,3] 

DigglePlots(Matern.pts, Matern.times)
paste("Matern II",toString(lambda),gsub("\\.","p",toString(r)),toString(nrow(Matern.pts)),sep="_")
paste("Realization of a Matern II point process. Intensity ", toString(lambda), ". Minimum distance ", toString(r), ". Number of points: ", toString(nrow(Matern.pts)), ".", sep="")
paste("Diagnostic plots for a realization of a Matern II point process.  Intensity ", toString(lambda), ". Minimum distance ", toString(r), ". Number of points: ", toString(nrow(Matern.pts)), ".", sep="")





# Baltimore data
balt <- readRDS("balt_standardized.rds")
balt.pts <- as.points(list(x = balt$coords[,1], y = balt$coords[,2]))
balt.times <- balt$times

plot(pp3(balt.pts[,1], balt.pts[,2], balt.times,as.box3(c(0,1,0,1,0,1))), 
     main = "Baltimore", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

DigglePlots(balt.pts, balt.times)


# Baltimore data with no time change
times.u <- unique(balt$times)
coords.u <- unique(balt$coords)

plot(pp3(coords.u[,1], coords.u[,2],rep(0.1,nrow(coords.u)), as.box3(c(0,1,0,1,0,1)) ), 
     main = "Baltimore (all)", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

coords.all <- matrix(rep(c(t(coords.u)), length(times.u)), ncol=2, byrow = TRUE)
times.all <- rep(unique(balt$times),each = nrow(coords.u))

plot(pp3(coords.all[,1], coords.all[,2], times.all, as.box3(c(0,1,0,1,0,1))),
     main = "Baltimore (all)", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

DigglePlots(coords.all, times.all)


## Questions
# Why is permutation test valid
# Assumption of stationarity






## Limitations for Baltimore dataset
# time is hardly continuous
# nonstationarity?
# leaving out zero observations altogether 