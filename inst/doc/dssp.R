## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(sp)
library(gstat)
library(DSSP)
library(ggplot2)

data("meuse.all")
data("meuse")

## ----eval=TRUE,include=TRUE---------------------------------------------------
meuse.train <- meuse.all[1:155, ]
meuse.valid <- meuse.all[156:164, ]
coordinates(meuse.train) <- ~ x + y
coordinates(meuse.valid) <- ~ x + y

## ----eval=TRUE,include=TRUE---------------------------------------------------
N <- 10000 ## number of samples to draw from the DSSP model
meuse.fit <- DSSP(
  formula = log(zinc) ~ 1, data = meuse.train, N = N, 
  pars = c(0.001, 0.001), log_prior = function(x) -2 * log(1 + x)
)

ETA <- meuse.fit$eta
DELTA <- meuse.fit$delta
NU <- meuse.fit$nu

## ----eval=TRUE,include=TRUE,fig.align="center",fig.width=5.5,fig.height=4.25----
Yhat <- rowMeans(exp(meuse.fit$y_fitted))

meuse$Yhat <- Yhat ## Model estimates of E(zinc concentration (ppm))
meuse$Y.true <- meuse.all$zinc[1:155]

##  Compare the smoothed values and the observed values

smooth.data <- data.frame(Yhat = meuse$Yhat, Y.true = meuse$Y.true)

smooth.scatterplot <- ggplot(smooth.data, aes(x = Yhat, y = Y.true)) +
  geom_point(size = 3) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Smoothed Values", y = "Observed Values", title = "Smoothed vs. Observed Values") +
  xlim(min(smooth.data), max(smooth.data)) +
  ylim(min(smooth.data), max(smooth.data)) +
  theme(plot.title = element_text(hjust = 0.5))

smooth.scatterplot

## ---- eval=TRUE,include=TRUE--------------------------------------------------
##  Now plot Parameter Estimates and ACF Plots

eta.densityplot <- ggplot(data.frame(x = ETA)) +
  geom_density(aes(x = x)) +
  labs(x = expression(eta), y = "posterior density", title = expression("Posterior Density of " * eta)) +
  theme(plot.title = element_text(hjust = 0.5))

delta.densityplot <- ggplot(data.frame(x = DELTA)) +
  geom_density(aes(x = x)) +
  labs(x = expression(delta), y = "posterior density", title = expression("Posterior Density of " * delta)) +
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=TRUE,include=TRUE, fig.height=3, fig.width=3.475, fig.show='hold'----

eta.densityplot
delta.densityplot

## ---- eval=TRUE,echo=FALSE,include=TRUE---------------------------------------

eta_acf <- acf(ETA, plot = FALSE)
eta_acfdf <- with(eta_acf, data.frame(lag, acf))

eta.acfplot <- ggplot(data = eta_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(x = "Lag", y = "ACF", title = expression("ACF for Samples from Posterior of " * eta)) +
  theme(plot.title = element_text(hjust = 0.5))

delta_acf <- acf(DELTA, plot = FALSE)
delta_acfdf <- with(delta_acf, data.frame(lag, acf))

delta.acfplot <- ggplot(data = delta_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(x = "Lag", y = "ACF", title = expression("ACF for Samples from Posterior of " * delta)) +
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=TRUE,include=TRUE, fig.height=3, fig.width=3.475, fig.show='hold'----
eta.acfplot
delta.acfplot

## ----eval=TRUE,include=TRUE---------------------------------------------------
eta.cumsumplot <- ggplot(data.frame(x = 1:length(ETA), y = cumsum(ETA) / (1:length(ETA)))) +
  geom_line(aes(x = x, y = y)) +
  labs(x = "sample", y = expression(eta), title = bquote(atop("Cumuulative Mean of Samples", "from Posterior of" ~ eta))) +
  theme(plot.title = element_text(hjust = 0.5))



delta.cumsumplot <- ggplot(data.frame(x = 1:length(DELTA), y = cumsum(DELTA) / (1:length(DELTA)))) +
  geom_line(aes(x = x, y = y)) +
  labs(x = "sample", y = expression(eta), title = bquote(atop("Cumuulative Mean of Samples", "from Posterior of" ~ delta))) +
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=TRUE,include=TRUE, fig.height=3, fig.width=3.475, fig.show='hold'----
eta.cumsumplot
delta.cumsumplot

## ----eval=TRUE,include=TRUE---------------------------------------------------
Y.pred <- predict(meuse.fit, meuse.valid)
Y.pred <- exp(Y.pred)

pred.data <- data.frame(Yhat.pred = rowMeans(Y.pred), Y.true = meuse.all$zinc[156:164])

pred.scatterplot <- ggplot(pred.data, aes(x = Yhat.pred, y = Y.true)) +
  geom_point(size = 3) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "Predicted Values", y = "True Values", title = "Predicted vs. True Values") +
  xlim(min(pred.data), max(pred.data)) +
  ylim(min(pred.data), max(pred.data)) +
  theme(plot.title = element_text(hjust = 0.5))

pred.boxplot <- ggplot(stack(as.data.frame(t(Y.pred)))) +
  geom_boxplot(aes(x = ind, y = values)) +
  geom_point(data = data.frame(Y.true = meuse.all$zinc[156:164]), aes(x = 1:9, y = Y.true), shape = 4, size = 3) +
  labs(x = "", y = "Y", title = bquote(atop("Boxplot of Predicted Values of", "Y and True Values (X)"))) +
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=TRUE,include=TRUE, fig.height=3, fig.width=3.475, fig.show='hold'----
pred.scatterplot
pred.boxplot

