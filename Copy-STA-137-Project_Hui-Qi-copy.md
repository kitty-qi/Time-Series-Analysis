Project
================
Hui Qi
3/13/2022

## R Markdown

``` r
library(readxl)
TempNH <- read_excel("TempNH_1850_2021.xlsx", col_names = c('x','y'))
x = TempNH[,1]
y = TempNH[,2]
y = unlist(y)
```

``` r
summary(y)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -0.70200 -0.36125 -0.12400 -0.03456  0.10075  1.27600

``` r
var(y)
```

    ## [1] 0.1901463

``` r
tm = 1:172
loesstrnd=loess(y~tm, span  = 0.25)

plot(tm, y, type="l", lty=1, xlab="Time", ylab="Annual Temperature Anomalies ", main="Time series with loess")
points(tm, loesstrnd$fitted, type="l", lty=2, col= "red")
legend("topleft", "loess", lty = 2, col = "red")
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
rough = loesstrnd$residuals
plot(rough, type='l',xlab = 'Time', ylab='Residuals', main = 'Rough part')
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
hist(loesstrnd$residuals, freq = F, xlab="residuals", main="Loess:histogram of residuals")
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
sse <- sum((loesstrnd$residuals)^2)
ssto <- sum((y-mean(y))^2)
R2 <- 1-sse/ssto
R2
```

    ## [1] 0.9264963

``` r
qqnorm(rough, main = "Normal probability plot of Rough")
qqline(rough)
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
par(mfrow =c(2,2))
plot.ts(y, main = 'Plot of the time series')

y1 = diff(y,1)
plot.ts(y1, main = 'First diff of ')

par(mfrow =c(1,1))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
par(mfrow =c(2,2))
acf(y1, main = 'ACF plot of first diff')
pacf(y1, main = 'PACF plot of first diff')

par(mfrow =c(1,1))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
library(astsa)
model1 = sarima(y, p=3,d=1,q=1,details=FALSE)

model1res= model1$fit$residuals

plot(model1res, type='l',xlab = 'Time', ylab='Residuals', main = 'Residuals of ARIMA(3,1,1)')
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
par(mfrow=c(2,2))
hist(model1res, freq = F, xlab="residuals", main="Histogram of residuals of ARIMA(3,1,1)")

qqnorm(model1res, main = "Normal probability plot of Residuals")
qqline(model1res)

par(mfrow=c(1,1))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
AIC = matrix(0,4,4)
for (i in 1:4){
  for (j in 1:4){
    AIC[i,j] <- sarima(y, p = i-1, d=1, q = j-1, details=FALSE)$AIC
  }
}
```

    ## Warning in sqrt(diag(fitit$var.coef)): NaNs produced

    ## Warning in sqrt(diag(fitit$var.coef)): NaNs produced

``` r
AIC
```

    ##            [,1]      [,2]      [,3]      [,4]
    ## [1,] -0.9229414 -1.094005 -1.119878 -1.108348
    ## [2,] -0.9922384 -1.114102 -1.108236 -1.109690
    ## [3,] -1.0455561 -1.114930 -1.112556 -1.114580
    ## [4,] -1.1233958 -1.116518 -1.110752 -1.103335

``` r
# -1.1233958 is the smallest ARIMA(3,1,0)

model2<-arima(y,order=c(3,1,0))
acf(model2$residuals)
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
par(mfrow=c(2,2))
hist(model2$residuals, freq = F, xlab="residuals", main="Histogram of residuals of ARIMA(3,1,0)")

qqnorm(model2$residuals, main = "Normal probability plot of Residuals")
qqline(model1res)

par(mfrow=c(1,1))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
model2
```

    ## 
    ## Call:
    ## arima(x = y, order = c(3, 1, 0))
    ## 
    ## Coefficients:
    ##           ar1      ar2      ar3
    ##       -0.4113  -0.3430  -0.2837
    ## s.e.   0.0738   0.0759   0.0739
    ## 
    ## sigma^2 estimated as 0.01821:  log likelihood = 99.61,  aic = -191.21

``` r
specselect=function(y,kmax){
# Obtains the values of the criterion function for
# obtaining the optimal number of neighbors for
# spectral density estimate for modified Daniell's method.
# input: y, observed series; kmax=max number of neighbors to
# be considered
# output: ctr - the criterion function
# output: kopt - the value of k at which the criterion function # is minimized
ii=spec.pgram(y,log="no",plot=FALSE)
ii=ii$spec
cc=norm(as.matrix(ii),type="F")^2 
ctr=rep(1,kmax) ###criterion function

for(k in 1:kmax) {
ss=2*k+1; kk=1/(2*k) 
ff=spec.pgram(y,spans=ss,log="no",plot=FALSE) 
fspec=ff$spec 
ctr[k]=norm(as.matrix(ii-fspec),type="F")^2+kk*cc
}

kopt=which.min(ctr)
result=list(ctr=ctr,kopt=kopt)
return(result)
}

specselect(y1,12)
```

    ## $ctr
    ##  [1] 0.07632212 0.06450387 0.06076200 0.05949365 0.05867776 0.05712129
    ##  [7] 0.05707398 0.05688819 0.05508067 0.05457552 0.05466916 0.05429911
    ## 
    ## $kopt
    ## [1] 12

``` r
plot(c(1:12),specselect(y1,12)$ctr,type="o")
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
koptimal<-specselect(y1,12)$kopt ##the one which minimizes the criterion function 
spans<-koptimal*2+1 ##optimal span
spans
```

    ## [1] 25

``` r
model2
```

    ## 
    ## Call:
    ## arima(x = y, order = c(3, 1, 0))
    ## 
    ## Coefficients:
    ##           ar1      ar2      ar3
    ##       -0.4113  -0.3430  -0.2837
    ## s.e.   0.0738   0.0759   0.0739
    ## 
    ## sigma^2 estimated as 0.01821:  log likelihood = 99.61,  aic = -191.21

``` r
spec_smooth <- spec.pgram(y1, spans=25, log="no", plot=FALSE)
freq <- spec_smooth$freq
spec <- spec_smooth$spec

library(astsa)
arma.spec(ar = model2$coef[1:3], var.noise = model2$sigma2, ylim = c(0,0.07))
lines(freq,spec, col = "blue")
legend("top", legend = c("spectral density function", "smoothed periodogram"), lty = c(1,1), col = c("black", "blue"))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
n = length(y)
ynew<- y[1:(n-6)]
ylast<- y[(n-5):n]

model4<-arima(ynew, order = c(3,1,0))
h <- 6
m <- n-h
fcast <- predict(model4, n.ahead=h)
upper <- fcast$pred+1.96*fcast$se
lower <- fcast$pred-1.96*fcast$se
#plot
plot.ts(ynew, xlim = c(0,n), xlab = "y", ylim=c(-0.7,1.4))
polygon(x=c(m+1:h,m+h:1), y=c(upper, rev(lower)), col="lightblue",border=NA)
lines(x=m+(1:h), y=fcast$pred, col="blue")
lines(x=m+(1:h), y=ylast, col="black")
legend("top", legend = c("true","fitted"), lty=c(1, 1), col = c("black","blue"))
```

![](Copy-STA-137-Project_Hui-Qi-copy_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
