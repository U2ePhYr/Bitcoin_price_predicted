library(readr)
library(TSA)
library(fUnitRoots)
library(forecast)
library(CombMSC)
library(lmtest)
library(fGarch)
library(rugarch)
library(AID)
library(zoo)
library(dplyr)
library(nortest)
library(FSAdata)
library(tseries)
library(readxl)



bitcoin <- read.csv("~/desktop/Bitcoin_Historical_Price.csv", header = TRUE)
days10price <- read_excel("~/desktop/Bitcoin_Prices_Forecasts.xlsx")

bitcoin$Close <- as.numeric(as.character(gsub(",","",bitcoin$Close)))
class(bitcoin)
head(bitcoin)
class(days10price)

# Convert data into a time series object
bitcoin.ts = matrix(bitcoin$Close, nrow = 2130, ncol = 1)
bitcoin.ts = as.vector(t(bitcoin.ts))
bitcoin.ts = ts(bitcoin$Close,frequency=365, start=c(2013,117), end = c(2019,55))
class(bitcoin.ts)

days10.ts = matrix(days10price$`Closing price`, nrow = 10, ncol = 1)
days10.ts = as.vector(t(days10.ts))
days10.ts = ts(days10price$`Closing price`,frequency=365, start=c(2019,56), end = c(2019,65))
class(days10.ts)

raw <- bitcoin.ts
raw1 <- days10.ts

sum(is.na(bitcoin.ts)) # Identify the missing value
# No Na in this dataset

# Load some useful functions:
MASE = function(observed , fitted ){
  # observed: Observed series on the forecast period
  # fitted: Forecast values by your model
  Y.t = observed
  n = length(fitted)
  e.t = Y.t - fitted
  sum = 0 
  for (i in 2:n){
    sum = sum + abs(Y.t[i] - Y.t[i-1] )
  }
  q.t = e.t / (sum/(n-1))
  MASE = data.frame( MASE = mean(abs(q.t)))
  return(list(MASE = MASE))
}

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

BoxCoxSearch = function(y, lambda=seq(-3,3,0.01), 
                        m= c("sf", "sw","ad" ,"cvm", "pt", "lt", "jb"), plotit = T, verbose = T){
  N = length(m)
  BC.y = array(NA,N)
  BC.lam = array(NA,N)
  for (i in 1:N){
    if (m[i] == "sf"){
      wrt = "Shapiro-Francia Test"
    } else if (m[i] == "sw"){
      wrt = "Shapiro-Wilk  Test"
    } else if (m[i] == "ad"){
      wrt = "Anderson-Darling Test"
    } else if (m[i] == "cvm"){
      wrt = "Cramer-von Mises Test"
    } else if (m[i] == "pt"){
      wrt = "Pearson Chi-square Test"
    } else if (m[i] == "lt"){
      wrt = "Lilliefors Test"
    } else if (m[i] == "jb"){
      wrt = "Jarque-Bera Test"
    } 
    
    print(paste0("------------- ",wrt," -------------"))
    out = tryCatch({boxcoxnc(y, method = m[i], lam = lambda, lambda2 = NULL, plot = plotit, alpha = 0.05, verbose = verbose)
      BC.lam[i] = as.numeric(out$lambda.hat)}, 
      error = function(e) print("No results for this test!"))
    
  }
  return(list(lambda = BC.lam,p.value = BC.y))
}

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH")[1]){
  # If you have an output from arima() function use class = "ARIMA"
  # If you have an output from garch() function use class = "GARCH"
  # If you have an output from ugarchfit() function use class = "ARMA-GARCH"
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  acf(res.model,main="ACF of standardised residuals")
  pacf(res.model,main="PACF of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
}

plot(bitcoin.ts,ylab='Bitcoin Price',xlab='Date',type='o', main = 'Time series plot of Bitcoin price daily')
# There are Changing variance and trend obviously. Meanwhile, intervention and 
# autoregressive behavior should be here. Furthermore, seasonality is not obviously 
# but fluctuating behavior probably.

par(mfrow=c(1,2))
acf(bitcoin.ts)
pacf(bitcoin.ts)
# Slowly decayin pattern in ACF and very high first correlation in PACF 
# implies the existence of trend and nonstationarity. Moreover, there are many
# high autocorrealation in PACF plot. Which implies there shoud exist changing 
# variance but no seasonality.
# ACF, PACF, EACF all shows the ARCH effect in this series.

bc.search = BoxCoxSearch(y=bitcoin.ts, lam=seq(-3,3,0.01), m= c("sf", "sw","ad" ,"cvm", "pt", "lt", "jb","ac"), plotit = T, verbose = T)
lambda = 0.30 # Applied for Box-Cox transformation for bitcoin series
BC.bitcoin.ts = ((bitcoin.ts^lambda)-1)/lambda

par(mfrow=c(1,1))
qqnorm(BC.bitcoin.ts)
qqline(BC.bitcoin.ts, col = 2)
shapiro.test(BC.bitcoin.ts)
# The Box-Cox Transformation did not help to improve normality of the series 
# obviously. Because many dots are not aligned with the red line in QQ plot 
# and p-value of the Shapiro test is less than 0.05.

ar(diff(BC.bitcoin.ts))
adfTest(BC.bitcoin.ts, lags = 33,  title = NULL,description = NULL)
# Adftest shows that this series is non-staionary after applied for Box-Cox 
# transformation.

diff.bitcoin.ts.BC = diff(BC.bitcoin.ts,differences = 1) # Do the first differencing for it
ar(diff(diff.bitcoin.ts.BC))
adfTest(diff.bitcoin.ts.BC, lags = 33,  title = NULL,description = NULL)
# The results of adfTest implies that the series is staionary after applying 
# the first differencing for Box-Cox transformed series.

par(mfrow=c(1,1))
plot(diff.bitcoin.ts.BC,type='o',main='Time series plot of Bitcoin price daily for the first differencing',ylab='Bitcoin Price')
# Time series plot, ACF and PACF do not show a clear trend. So I'll go on with the first difference.

par(mfrow=c(1,2))
acf(diff.bitcoin.ts.BC)
pacf(diff.bitcoin.ts.BC)
# Both ACF and PACF shows that more than 5 high correlation exist separately.
# What's more, it seems to a pattern in ACF plot.
# So i can consider ARIMA(5,1,0) and ARIMA(5,1,5) for possible models.

eacf(diff.bitcoin.ts.BC,ar.max = 10, ma.max = 10)
# Almost the whole fourth column in EACF are filled by the value of x, which
# means the ARCH effect is still in this series.
# From the output of eacf we include ARIMA(0,1,1), ARIMA(1,1,1), ARIMA(1,1,0)
# models in the set of possible models.
# ACF, PACF, EACF all shows the ARCH effect in this series.

par(mfrow=c(1,1))
res = armasubsets(y=diff.bitcoin.ts.BC,nar=5,nma=5,y.name='test',ar.method='ols')
plot(res)
# In the BIC table, the final set of possible models is
# {ARIMA(4,1,4), ARIMA(5,1,4), ARIMA(1,1,4)}

# Model Diagnostics: Because of non-normality for this BC series, only CSS 
# method be used below.

# ARIMA(0,1,1)
model_011_css = arima(BC.bitcoin.ts,order=c(0,1,1),method='CSS')
coeftest(model_011_css)
# All coefficients are significant. This could be considered for good model.

# ARIMA(1,1,1)
model_111_css = arima(BC.bitcoin.ts,order=c(1,1,1),method='CSS')
coeftest(model_111_css)
# Both coefficients are insignificant. So it should be removed from possible 
# models.

# ARIMA(1,1,0)
model_110_css = arima(BC.bitcoin.ts,order=c(1,1,0),method='CSS')
coeftest(model_110_css)
# All coefficients are insignificant. So it should be removed from possible 
# models.

# ARIMA(4,1,4)
model_414_css = arima(BC.bitcoin.ts,order=c(4,1,4),method='CSS')
coeftest(model_414_css)
# All coefficients are significant. This could be considered for good model.

# ARIMA(5,1,4)
model_514_css = arima(BC.bitcoin.ts,order=c(5,1,4),method='CSS')
coeftest(model_514_css)
# Most of the coefficients are insignificant. So it should be removed from possible 
# models.

# ARIMA(1,1,4)
model_114_css = arima(BC.bitcoin.ts,order=c(1,1,4),method='CSS')
coeftest(model_114_css)
# Most of the coefficients are insignificant. So it should be removed from possible 
# models.

# ARIMA(5,1,0)
model_510_css = arima(BC.bitcoin.ts,order=c(5,1,0),method='CSS')
coeftest(model_510_css)
# For the same reason, i will reremove it from possible models.

# ARIMA(5,1,5)
model_515_css = arima(BC.bitcoin.ts,order=c(5,1,5),method='CSS')
coeftest(model_515_css)
# For the same reason, i will reremove it from possible models.

# Now i got two possible good models ARIMA(0,1,1) and ARIMA(4,1,4), then do 
# overdifferencing Diagnostics for these two model.

# ARIMA(0,1,2)
model_012_css = arima(BC.bitcoin.ts,order=c(0,1,2),method='CSS')
coeftest(model_012_css)
# Both coefficients are insignificant. So it should be removed from possible 
# models.

# ARIMA(4,1,5)
model_415_css = arima(BC.bitcoin.ts,order=c(4,1,5),method='CSS')
coeftest(model_415_css)
# Half of the coefficients are insignificant. So it should be removed from possible 
# models.

# Now i got two good models ARIMA(0,1,1) and ARIMA(4,1,4).

# Because of non-normality for this Box-Cox bitcoin series, only CSS 
# method be used. So i cant give BIC and AIC table for these two models.

graphics.off() 
residual.analysis(model_011_css, std = TRUE,start = 1) # Do residual analysis
# From the time series plot of standardised residuals, should be conclude 
# that there is still changing variance. Correspondingly, both results in 
# QQ plot which is two tails are represented and Ljung-Box Test are confirmed 
# ARCH effect again. Furthermore, both ACF and PACF imply that there should 
# be correlation problem in series. One good thing is symmetric histogram 
# been get.   

graphics.off() 
residual.analysis(model_414_css, std = TRUE,start = 1)
# A little bit better than what's the performance in model_011's residual analysis.
# Especially for autocorrelation parts. However, both results of QQ plot 
# and Ljung-Box Test indicate ARCH effect in this series.

# Go for fit GARCH model.
res.m.414 = residuals(model_414_css)

par(mfrow=c(1,1))
plot(res.m.414,type='o',ylab="Bitcoin Price",main="Time series plot of Bitcoin price residuals daily")
par(mfrow=c(1,2))
acf(res.m.414)
pacf(res.m.414)
eacf(res.m.414)
# From the model ARIMA414 residual series of bitcoin price, there is an ARCH effect present 
# in the series while it is staionary and its conditional mean is zero. 
# Which represents by volatility clustering behavior. Meanwhile,
# Both ACF and PACF shows the ARCH effect in this series.

# Generally this series should be considered about fitting ARIMA+GARCH models

abs.res.m.414 = abs(res.m.414)
sqr.res.m.414 = res.m.414^2

par(mfrow=c(1,2))
acf(abs.res.m.414)
pacf(abs.res.m.414)
eacf(abs.res.m.414,ar.max = 5, ma.max = 5)
# Both ACF and PACF show slowly decayin pattern. So got from eacf  
# {ARMA(2,2), ARMA(1,2), ARMA(1,3)} ==> {GARCH(2,2), GARCH(2,1), GARCH(3,1)}

par(mfrow=c(1,2))
acf(sqr.res.m.414)
pacf(sqr.res.m.414)
eacf(sqr.res.m.414,ar.max = 5, ma.max = 5)
# Both ACF and PACF show slowly decayin pattern. So got from eacf
#{ARMA(1,2), ARMA(1,3)} ==> {GARCH(2,1), GARCH(3,1)}

model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), 
                   mean.model = list(armaOrder = c(0,1,1), include.mean = FALSE), 
                   distribution.model = "norm")
m1<-ugarchfit(spec=model1,data=res.m.414)
m1  # AIC = 0.3496
plot(m1,which=2)
# Only two coefficients are insignificant. So this is not bad model in terms of 
# significant test. From the plot of model1, the Density of standardized 
# residuals looks good. Meanwhile only ACF of squared standardized residuals 
# is acceptable. However, the assumption of normality is not good.

model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 1)), 
                   mean.model = list(armaOrder = c(0,1,1), include.mean = FALSE), 
                   distribution.model = "norm")
m2<-ugarchfit(spec=model2,data=res.m.414)
m2  # AIC = 0.3589
plot(m2,which=2)
# The result of significant test is similar with model1. From the plot of model2, 
# the Density of standardized residuals looks good. Meanwhile only ACF of 
# squared standardized residuals is acceptable. However, the assumption of 
# normality is not good.

model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(3, 1)), 
                   mean.model = list(armaOrder = c(0,1,1), include.mean = FALSE), 
                   distribution.model = "norm")
m3<-ugarchfit(spec=model3,data=res.m.414)
m3  # AIC = 0.3599
plot(m3,which=2)
#...

model4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), 
                   mean.model = list(armaOrder = c(4,1,4), include.mean = FALSE), 
                   distribution.model = "norm")
m4<-ugarchfit(spec=model4,data=res.m.414)
m4  # AIC = 0.3447(The least)
plot(m4,which=2)
# Only three coefficients are insignificant in totally ten coefficients. 
# So this is a good model in terms of significant test. From the plot of model1, 
# the Density of standardized residuals looks good. Meanwhile both ACF of 
# squared standardized residuals and standardized residuals are acceptable almostly.
# However, the assumption of normality is not good.

model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 1)), 
                   mean.model = list(armaOrder = c(4,1,4), include.mean = FALSE), 
                   distribution.model = "norm")
m5<-ugarchfit(spec=model5,data=res.m.414)
m5  # AIC = 0.3549
plot(m5,which=2)
#....

model6<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(3, 1)), 
                   mean.model = list(armaOrder = c(4,1,4), include.mean = FALSE), 
                   distribution.model = "norm")
m6<-ugarchfit(spec=model6,data=res.m.414)
m6  # AIC = 0.3558
plot(m6,which=2)
#....
# Model 4 is the best model above in terms of AIC.
model_414_css = arima(bitcoin.ts, order=c(4,1,4),method='CSS')
res.414 = residuals(model_414_css)

m.22 = garch(res.414,order=c(2,2),trace = FALSE)
m.22_2 = garchFit(formula = ~garch(2,2), data =res.414, algorithm = "lbfgsb" )
summary(m.22_2)
res = m.22_2@residuals
par(mfrow=c(1,2))
acf(res)
pacf(res)
par(mfrow=c(1,1))
plot(m.22_2@residuals)
# From the results of Ljung-Box Test, the data are independently distributed.
# Which means there is no correlation in residuals datasets.
# Meanwhile, both ACF and PACF plots are white noise almostly.

# For fitted valuea:
m.22_2 <- ugarchspec() # Here I assume your model object is called "model"
m.fit <-ugarchfit(spec=m.22_2,data=raw)
fitted.values = fitted(m.fit)

# MASE for fitted value:
MASE(raw, fitted.values)

# For forecast values:
m.22_2<-ugarchspec()
m.fit <-ugarchfit(spec=m.22_2,data=raw)
forc = ugarchforecast(m.fit, data = raw, n.ahead = 10)
forecasts = forc@forecast$seriesFor

# MASE for forecast value:
MASE(raw1, forecasts)
































