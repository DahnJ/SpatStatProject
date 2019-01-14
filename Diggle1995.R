### Analysis of Diggle's 1995 paper ###

library(splancs)
library(spatstat)


# Function for automatic plots. See ?stmctest for explanation of parameters
# Returns the output of stmctest
DigglePlots <- function(pts, times, poly, tlimits = c(0,1), s, tt, nsim = 100, quiet = FALSE, Dzero = FALSE){
  if(missing(poly)){ poly <- bboxx(spoints(c(0,1,0,1)))}
  if(missing(s)){ s <- seq(0,0.7,0.01)[-1]}
  if(missing(tt)) { tt <- s }
  
  stkh <- stkhat(pts, times, poly, tlimits, s, tt)
  stse <- stsecal(pts, times, poly, tlimits, s, tt)
  stmc <- stmctest(pts, times, poly, tlimits, s, tt, nsim, quiet)
  
  stdiagn(pts, stkh, stse, stmc, Dzero)
  return(stmc)
}



## Poisson
lambda <- 1000
poisson <- rpoispp3(lambda)
plot(poisson, main = "Poisson", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)
distances <- seq(0,0.7,0.01)[-1]
times <- distances
poisson.pts <- as.points(list(x = poisson$data$x, y = poisson$data$y))
poisson.times <- poisson$data$z

MC.Poisson <- DigglePlots(poisson.pts, poisson.times)

paste("Poisson",toString(lambda),toString(nrow(poisson.pts)),sep="_")
paste("Realization of a Poisson point process. Intensity ", str(lambda), ". Number of points: ", toString(nrow(poisson.pts)), sep="")
paste("Diagnostic plots for a realization of a Poisson point process. Intensity ", str(lambda), ". Number of points: ", toString(nrow(poisson.pts)), sep="")



## Multidim Thomas process
lambdap <- 20
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


lambdac <- 40
sigmac <- 0.02

Thomas <- apply(points,1,generateChildren, lambda = lambdac, sigma = sigmac)
Thomas <- do.call(rbind,Thomas)

Thomas.pts <- as.points(list(x = Thomas[,1], y = Thomas[,2]))
Thomas.times <- Thomas[,3]

plot(pp3(Thomas.pts[,1],Thomas.pts[,2],Thomas.times,as.box3(c(0,1,0,1,0,1))), 
     main = "Thomas", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

MC.Thomas <- DigglePlots(Thomas.pts, Thomas.times)

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

MC.Mattern <- DigglePlots(Matern.pts, Matern.times)

paste("Matern II",toString(lambda),gsub("\\.","p",toString(r)),toString(nrow(Matern.pts)),sep="_")
paste("Realization of a Matern II point process. Intensity ", toString(lambda), ". Minimum distance ", toString(r), ". Number of points: ", toString(nrow(Matern.pts)), ".", sep="")
paste("Diagnostic plots for a realization of a Matern II point process.  Intensity ", toString(lambda), ". Minimum distance ", toString(r), ". Number of points: ", toString(nrow(Matern.pts)), ".", sep="")





## Baltimore data
balt <- readRDS("balt_standardized.rds")
balt.pts <- as.points(list(x = balt$coords[,1], y = balt$coords[,2]))
balt.times <- balt$times

plot(pp3(balt.pts[,1], balt.pts[,2], balt.times,as.box3(c(0,1,0,1,0,1))), 
     main = "Baltimore", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

MC.balt <- DigglePlots(balt.pts, balt.times)


# Baltimore data with no time change
times.u <- unique(balt$times)
coords.u <- unique(balt$coords)

# plot(pp3(coords.u[,1], coords.u[,2],rep(0.1,nrow(coords.u)), as.box3(c(0,1,0,1,0,1)) ), 
#     main = "Baltimore (all)", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

coords.all <- matrix(rep(c(t(coords.u)), length(times.u)), ncol=2, byrow = TRUE)
times.all <- rep(unique(balt$times),each = nrow(coords.u))

plot(pp3(coords.all[,1], coords.all[,2], times.all, as.box3(c(0,1,0,1,0,1))),
     main = "Baltimore (all)", box.front=list(lty=5,col="red"), box.back=list(lty=5,col="red"), cex=0.8)

MC.baltall <- DigglePlots(coords.all, times.all)