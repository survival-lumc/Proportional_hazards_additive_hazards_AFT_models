#' ---
#' title: "An illustration of proportional hazards, additive hazards
#' and accelerated failure time models"
#' author: "Hein Putter"
#' date: "`r Sys.setenv(LANG = 'en_US.UTF-8'); format(Sys.Date(), '%d %B %Y')`"
#' output:
#'   pdf_document:
#'     toc: true
#'   latex_engine: xelatex
#' fontsize: 12pt
#' ---
#' 
#' # Preliminaries
#' 
#' This document serves as the Supplementary Material of the
#' paper "Alternatives to the Cox proportional hazards model: A
#' review of additive hazards and accelerated failure time
#' models", (to be) published in Statistica Neerlandica. Section
#' 5.2 of that paper contains a data example, of which this
#' document shows the full R code.
#' 
#' We will be using the pbc data from the survival package, in a
#' very much simplified setting. The aim of this illustration is
#' to highlight differences in assumptions between different
#' models, and to illustrate techniques that attempt to alleviate
#' restrictions in the original models. We therefore restrict to
#' only two covariates, hepato and edema, out of many in the pbc
#' data that significantly influence survival. Of these two,
#' edema violates the PH assumption.
# Suppresses all warnings for the rest of this script (not just
# the package-loading messages below) - convenient here since
# several packages used produce noisy convergence/deprecation
# warnings, but it also hides any other warning that might occur.
options(warn = -1)
# Loading packages
suppressPackageStartupMessages({
  library(tictoc) # used for computation time
  library(survival)
  library(ggplot2) # used for all survival-curve plots
  library(patchwork) # used for the combined figure at the end
  # Packages for additive hazards models, in addition to survival
  library(ahMLE)
  library(timereg)
  library(aftgee)
  # Packages for AFT models, in addition to survival
  library(lss2)
  library(rms)
  library(aftsem)
  library(flexsurv)
  library(eha)
})

# head(pbc) # look at help(pbc) for more information on the data
table(pbc$hepato, useNA = "always")
#' One of the two variables of interest, `hepato`, is only used
#' in the original 312 patients included in the RCT (see
#' `help(pbc)` for context). For the purpose of this illustration
#' we restrict to this subset. We are going to overwrite pbc (the
#' original pbc data will remain available in the survival
#' package).
pbc <- subset(pbc, !is.na(hepato))
# Looking at events
table(pbc$status) # 0=censored, 1=transplant, 2=dead
pbc$stat <- pbc$status
pbc$stat[pbc$stat == 2] <- 1 # so we use transplant-free survival
table(pbc$stat)
table(pbc$edema)
#' We are also going to overwrite edema with a binary version
#' where the value 0.5 is set to 1 (so 0.5 and 1 are combined
#' into one group). (Again, the original data will remain
#' available in the survival package.)
pbc$edema[pbc$edema == 0.5] <- 1
#' In what follows we are going to construct model-based survival
#' curves for each of the four combinations of `hepato` and
#' `edema`. These are the numbers for these combinations.
table(pbc$hepato, pbc$edema)
# We prefer to work with time in years
year <- 365.25
pbc$yrs <- pbc$time / year

# Helper for the four-group (hepato x edema) survival curve plots
# used throughout this script in place of base R plot()/lines()/
# legend(). type = "step" for Kaplan-Meier-type/model-based
# curves defined by jumps at event times; type = "line" for
# smooth parametric curves evaluated on a dense time grid. No
# confidence bands or numbers-at-risk are shown, per the level of
# detail used in the original base R plots. combine_he() and
# plot_he_survival() below both rely on the global constants
# he_levels and he_pal defined here.
he_levels <- c("No / no", "No / yes", "Yes / no", "Yes / yes")
# matches base R col = 1:4
he_pal <- c("black", "red", "green3", "blue")

# Combines the four per-group (hepato x edema) data frames - each
# with columns time and surv - into one long data frame with a
# group column, ready for plot_he_survival().
combine_he <- function(d00, d01, d10, d11) {
  d00$group <- he_levels[1]
  d01$group <- he_levels[2]
  d10$group <- he_levels[3]
  d11$group <- he_levels[4]
  out <- rbind(d00, d01, d10, d11)
  out$group <- factor(out$group, levels = he_levels)
  out
}

plot_he_survival <- function(
    df, title, type = c("step", "line")
) {
  type <- match.arg(type)
  geom_fun <- if (type == "step") geom_step else geom_line
  ggplot(df, aes(x = time, y = surv, color = group)) +
    geom_fun(linewidth = 0.8) +
    scale_color_manual(values = he_pal, breaks = he_levels) +
    # coord_cartesian() (rather than ylim()/scale limits) zooms
    # the view to [0, 1] without dropping or warning about data
    # points outside that range - needed because one of the
    # additive hazards curves later on goes slightly above 1.
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Years",
      y = "Survival",
      color = "Hepato / edema",
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

#' We start with separate Kaplan-Meier survival curves for the
#' four different combinations no/yes of hepato and edema.
pbc$he <- 2 * pbc$hepato + pbc$edema
pbc$he <- factor(
  pbc$he, levels = 0:3,
  labels = c("No / no", "No / yes", "Yes / no", "Yes / yes")
)
sfhe <- survfit(Surv(yrs, stat) ~ he, data = pbc)
# pdf("Kaplan_Meiers.pdf")
kmdf <- data.frame(
  time = sfhe$time,
  surv = sfhe$surv,
  group = rep(sub("^he=", "", names(sfhe$strata)), sfhe$strata)
)
kmdf$group <- factor(kmdf$group, levels = he_levels)
# plot.survfit() implicitly starts every curve at (0, 1); add
# that starting point explicitly since we are bypassing
# plot.survfit().
kmdf <- rbind(
  data.frame(time = 0, surv = 1,
             group = factor(he_levels, levels = he_levels)),
  kmdf
)
# pdf("Kaplan_Meiers.pdf")
plot_he_survival(
  kmdf, "Kaplan-Meier survival curves", type = "step"
)
# dev.off()

#' # Proportional hazards models
#' 
#' We start by fitting a Cox proportional hazards model with the
#' two variables of interest.
c1 <- coxph(Surv(yrs, stat) ~ hepato + edema, data = pbc)
summary(c1) # both highly significant
cz1 <- cox.zph(c1)
cz1 # edema violates PH
#' From the output of `cz1` we see that `edema` violates the
#' proportional hazards (PH) assumption. The time-dependent
#' behaviour $\beta(t)$ can be plotted, from which we can learn
#' how the estimated regression coefficients varies with time. I
#' am using `transform = "identity"` here to make sure that in
#' the plot we see time linearly. The default is to plot
#' $\beta(t)$ with time transformed, which makes it more
#' difficult (in my view) to interpret $\beta(t)$.
par(mfrow = c(1, 2))
plot(cox.zph(c1, transform = "identity"))
par(mfrow = c(1, 1))
#' The $\beta(t)$ for `hepato` seems to go up and down somewhat,
#' but one could imagine a straight line reasonably within the
#' 95\% pointwise confidence intervals, compatible with a
#' time-constant $\beta$. The $\beta(t)$ for `edema` is quite
#' different. The estimated $\beta$ is huge, about 4 (this is the
#' *log* hazard ratio!), and quickly decreases to 0 after about 5
#' years.
#' 
#' The easiest way to deal with the violation of the PH
#' assumption for edema is to fit a stratified Cox model, where
#' separate baseline hazards are estimated for edema no and edema
#' yes, and a hazard ratio (HR) for hepato (assumed to be the
#' same for edema no and edema yes).
# Stratified Cox model, stratified according to edema
c2 <- coxph(Surv(yrs, stat) ~ hepato + strata(edema), data = pbc)
summary(c2)
cox.zph(c2) # PH not violated for hepato
#' The PH assumption is not violated for hepato, so at least with
#' respect to the PH assumption this seems to be a reasonable
#' model. The model does not have an estimated hazard ratio /
#' regression coefficient for `edema`. The reason is that the PH
#' assumption for `edema` (which would give a HR) is replaced by
#' two separate baseline hazards (so no HR). If a numeric
#' description of the effect of `edema` is required, one would
#' have to model $\beta(t)$ by replacing the stratification by an
#' interaction of `edema` with a pre-specified function of time.
#' The plot of `cox.zph` above could inform which function to
#' use. Judging by this plot, something like $f(t) = \log(t + 1)$
#' could work well.
#' 
#' We will continue to show survival curves for four patients,
#' namely for all four combinations of hepato 0/1, and edema 0/1.
# According to Cox model, and stratified Cox model
# First Cox model
nd00 <- data.frame(hepato = 0, edema = 0)
nd01 <- data.frame(hepato = 0, edema = 1)
nd10 <- data.frame(hepato = 1, edema = 0)
nd11 <- data.frame(hepato = 1, edema = 1)

sf00 <- survfit(c1, newdata = nd00)
sf00 <- data.frame(time = sf00$time, surv = sf00$surv)
sf00 <- rbind(data.frame(time = 0, surv = 1), sf00)
sf00 <- subset(sf00, !duplicated(surv))
sf01 <- survfit(c1, newdata = nd01)
sf01 <- data.frame(time = sf01$time, surv = sf01$surv)
sf01 <- rbind(data.frame(time = 0, surv = 1), sf01)
sf01 <- subset(sf01, !duplicated(surv))
sf10 <- survfit(c1, newdata = nd10)
sf10 <- data.frame(time = sf10$time, surv = sf10$surv)
sf10 <- rbind(data.frame(time = 0, surv = 1), sf10)
sf10 <- subset(sf10, !duplicated(surv))
sf11 <- survfit(c1, newdata = nd11)
sf11 <- data.frame(time = sf11$time, surv = sf11$surv)
sf11 <- rbind(data.frame(time = 0, surv = 1), sf11)
sf11 <- subset(sf11, !duplicated(surv))

# pdf("Cox_model.pdf")
df_cox <- combine_he(sf00, sf01, sf10, sf11)
p_cox <- plot_he_survival(df_cox, "Cox model", type = "step")
p_cox
# dev.off()

sf00_cox <- sf00
sf01_cox <- sf01
sf10_cox <- sf10
sf11_cox <- sf11

# Then stratified Cox model
sf00 <- survfit(c2, newdata = nd00)
sf00 <- data.frame(time = sf00$time, surv = sf00$surv)
sf00 <- rbind(data.frame(time = 0, surv = 1), sf00)
sf00 <- subset(sf00, !duplicated(surv))
sf01 <- survfit(c2, newdata = nd01)
sf01 <- data.frame(time = sf01$time, surv = sf01$surv)
sf01 <- rbind(data.frame(time = 0, surv = 1), sf01)
sf01 <- subset(sf01, !duplicated(surv))
sf10 <- survfit(c2, newdata = nd10)
sf10 <- data.frame(time = sf10$time, surv = sf10$surv)
sf10 <- rbind(data.frame(time = 0, surv = 1), sf10)
sf10 <- subset(sf10, !duplicated(surv))
sf11 <- survfit(c2, newdata = nd11)
sf11 <- data.frame(time = sf11$time, surv = sf11$surv)
sf11 <- rbind(data.frame(time = 0, surv = 1), sf11)
sf11 <- subset(sf11, !duplicated(surv))

# pdf("Stratified_Cox_model.pdf")
df_scox <- combine_he(sf00, sf01, sf10, sf11)
p_scox <- plot_he_survival(
  df_scox, "Stratified Cox model", type = "step"
)
p_scox
# dev.off()

#' # Additive hazards models
#' 
#' We continue to fit additive hazards models. The non-parametric
#' OLS has been implemented for instance in the survival package,
#' function aareg(), and in the timereg package, using aalen().
#' We will first illustrate aareg() from survival.
a0 <- aareg(Surv(yrs, stat) ~ hepato + edema, data = pbc)
head(a0$coef) # Estimated beta's
head(apply(a0$coef, 2, cumsum)) # Estimated cumulative beta's
# plot(a0)
#' The function aalen() from the timereg package gives exactly
#' the same estimates as aareg() from survival, but in addition
#' useful tests for $\beta(t) = 0$ and for $\beta(t)$ being
#' time-constant are also provided (for the latter null
#' hypothesis both a Kolmogorov-Smirnov and a Cramer-von Mises
#' test).
a1 <- aalen(Surv(yrs, stat) ~ hepato + edema, data = pbc)
head(a1$cum) # same as the cumsum of a0$coef
summary(a1)
#' Both hepato and edema have a significant effect on survival
#' (null hypothesis of $\beta(t) = 0$ is rejected). The null
#' hypothesis of $\beta(t)$ being constant is rejected for edema,
#' not for hepato (note that this assumption is not the PH
#' assumption). A plot function for the result is provided.
# pdf("Ah_Bt.pdf", width = 8, height = 5)
par(mfrow = c(1, 3))
plot(a1, pointwise.ci = 2)
# dev.off()

#' We now derive the survival curves for the same four patients.
ah1 <- a1$cum
ah1 <- as.data.frame(ah1)
names(ah1) <- dimnames(ah1)[[2]]
names(ah1)[2] <- "intercept"
head(ah1)

b00 <- data.frame(time = ah1$time, Haz = ah1$intercept)
b00$surv <- exp(-b00$Haz)
b01 <- data.frame(
  time = ah1$time, Haz = ah1$intercept + ah1$edema
)
b01$surv <- exp(-b01$Haz)
b10 <- data.frame(
  time = ah1$time, Haz = ah1$intercept + ah1$hepato
)
b10$surv <- exp(-b10$Haz)
b11 <- data.frame(
  time = ah1$time,
  Haz = ah1$intercept + ah1$hepato + ah1$edema
)
b11$surv <- exp(-b11$Haz)

# pdf("Ah_nonpar_survival.pdf")
par(mfrow = c(1, 1))
df_ah1 <- combine_he(b00, b01, b10, b11)
p_ah1 <- plot_he_survival(
  df_ah1, "Additive hazards model (OLS)", type = "step"
)
p_ah1
# dev.off()

b00_OLS <- b00
b01_OLS <- b01
b10_OLS <- b10
b11_OLS <- b11

#' Note the non-monotone behavior of in particular the hepato no,
#' edema yes curve.
#' 
#' Since the effect of hepato seemed to be time-constant, we also
#' fit a model with a constant effect of hepato.
a2 <- aalen(Surv(yrs, stat) ~ const(hepato) + edema, data = pbc)
summary(a2)
par(mfrow = c(1, 2))
plot(a2, pointwise.ci = 2)

#' The estimated cumulative curves for intercept (cumulative
#' baseline hazard) and `edema` have not fundamentally changed,
#' but note that the cumulative baseline hazard estimate is
#' negative until just past $t=2$. Also for this model we derive
#' the survival curves for the same four patients.
ah2 <- a2$cum
ah2 <- as.data.frame(ah2)
names(ah2) <- dimnames(ah2)[[2]]
names(ah2)[2] <- "intercept"
head(ah2)

b00 <- data.frame(time = ah2$time, Haz = ah2$intercept)
b00$surv <- exp(-b00$Haz)
b01 <- data.frame(
  time = ah2$time, Haz = ah2$intercept + ah2$edema
)
b01$surv <- exp(-b01$Haz)
b10 <- data.frame(
  time = ah2$time,
  Haz = ah2$intercept + ah2$time * c(a2$gamma)
)
b10$surv <- exp(-b10$Haz)
b11 <- data.frame(
  time = ah2$time,
  Haz = ah2$intercept + ah2$time * c(a2$gamma) + ah2$edema
)
b11$surv <- exp(-b11$Haz)

# pdf("Ah_McKS_survival.pdf")
par(mfrow = c(1, 1))
df_ah2 <- combine_he(b00, b01, b10, b11)
p_ah2 <- plot_he_survival(
  df_ah2,
  "Additive hazards model (OLS)\ntime-constant hepato",
  type = "step"
)
p_ah2
# dev.off()

#' We can now also see the model-based survival curve for no
#' hepato, no edema being larger than one (consistent with the
#' cumulative baseline hazard being negative).

b00_OLS_tc <- b00
b01_OLS_tc <- b01
b10_OLS_tc <- b10
b11_OLS_tc <- b11

#' Finally we consider the MLE estimates of Lu et al. (2023).
a3 <- ah(Surv(yrs, stat) ~ hepato + edema, data = pbc)
a3 <- a3$cumbeta

# Survival curves for the same four patients
b00 <- data.frame(time = a3$time, Haz = a3$intercept)
b00$surv <- exp(-b00$Haz)
b01 <- data.frame(time = a3$time, Haz = a3$intercept + a3$edema)
b01$surv <- exp(-b01$Haz)
b10 <- data.frame(time = a3$time, Haz = a3$intercept + a3$hepato)
b10$surv <- exp(-b10$Haz)
b11 <- data.frame(
  time = a3$time,
  Haz = a3$intercept + a3$hepato + a3$edema
)
b11$surv <- exp(-b11$Haz)

# pdf("Ah_MLE_survival.pdf")
par(mfrow = c(1, 1))
df_ah3 <- combine_he(b00, b01, b10, b11)
p_ah3 <- plot_he_survival(
  df_ah3, "Additive hazards model (MLE)", type = "step"
)
p_ah3
# dev.off()

#' Note that all mode-based survival curves are monotone, but
#' especially the hepato no, edema yes curve is quite different
#' from the previous additive hazard curves.

b00_MLE <- b00
b01_MLE <- b01
b10_MLE <- b10
b11_MLE <- b11

#' # AFT models
#' 
#' Finally we consider AFT models, first parametric, then
#' semi-parametric AFT models.
#' 
#' ## Parametric AFT models
#' 
#' We start with parametric AFT models, first survreg() from the
#' survival package.
tic("survreg from {survival}")
ls.sur.weib <- survreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc
)
toc()
summary(ls.sur.weib)
#' It is instructive to see what happens if we use the original
#' time scale, days, rather than years. From the results we see
#' that the estimated regression coefficients for hepato and
#' edema remain the same, but the intercept has changed.
tic("survreg from {survival}, time in days")
ls.sur.weib2 <- survreg(
  Surv(time, stat) ~ hepato + edema, data = pbc
)
toc()
summary(ls.sur.weib2)
#' Many other parametric distributions can be used in survreg(),
#' like for instance the logistic distributions.
tic("survreg with logistic distribution")
ls.sur.logistic <- survreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "logistic"
)
toc()
summary(ls.sur.logistic)
#' About the parametrization of the Weibull distribution. The
#' function flexsurvreg() in the flexsurv package fits AFT
#' Weibull models with a parametrization more in line with that
#' used in the paper (and models with other parametric error
#' distributions). So let's fit the Weibull AFT model using
#' flexsurvreg(), and also AFT models with a range of other
#' distributions, including the generalized gamma distribution.
fitw <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "weibull"
)
fitl <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "llogis"
)
# Named fitgg (not fitg) to avoid colliding with the Gompertz fit
# below - the original script assigned both to fitg, so the
# gengamma fit was being silently overwritten by the Gompertz
# one.
fitgg <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "gengamma"
)
fitn <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "lnorm"
)
fitf <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "genf"
)
fitg <- flexsurvreg(
  Surv(yrs, stat) ~ hepato + edema, data = pbc,
  dist = "gompertz"
)
fitw
fitl
fitgg
fitn
fitf
fitg

#' Package eha has the function aftreg() that fits parametric AFT
#' models. Default is the Weibull distributions, but also
#' loglogistic, lognormal and Gompertz distributions are
#' supported. The same results are obtained as with other
#' packages, at least with respect to the regression
#' coefficients. It is unclear where for instance the shape and
#' scale of the Weibull model can be found.
aftregw <- aftreg(Surv(yrs, stat) ~ hepato + edema, data = pbc)
summary(aftregw)
aftregll <- aftreg(Surv(yrs, stat) ~ hepato + edema, data = pbc,
                   dist = "loglogistic")
summary(aftregll)

#' The function psm() from the rms package also fits AFT models.
#' From this it also seems to be easier to get predicted survival
#' functions. The following code provides model-based survival
#' curves for the four groups.
mySweibull <- function(x, rate, shape)
{
  return(exp( - rate * x^shape))
}
psm.weib <- psm(Surv(yrs, stat) ~ hepato + edema, data = pbc)
psm.weib
sig <- psm.weib$scale
lp <- predict(psm.weib, newdata = nd00, type = "lp")
tseq <- seq(0, 11, by = 0.01)
p00 <- mySweibull(tseq, shape = 1 / sig, rate = exp(-lp / sig))
lp <- predict(psm.weib, newdata = nd01, type = "lp")
p01 <- mySweibull(tseq, shape = 1 / sig, rate = exp(-lp / sig))
lp <- predict(psm.weib, newdata = nd10, type = "lp")
p10 <- mySweibull(tseq, shape = 1 / sig, rate = exp(-lp / sig))
lp <- predict(psm.weib, newdata = nd11, type = "lp")
p11 <- mySweibull(tseq, shape = 1 / sig, rate = exp(-lp / sig))

# pdf("AFT_Weibull_survival.pdf")
df_weib <- combine_he(
  data.frame(time = tseq, surv = p00),
  data.frame(time = tseq, surv = p01),
  data.frame(time = tseq, surv = p10),
  data.frame(time = tseq, surv = p11)
)
plot_he_survival(df_weib, "Weibull AFT", type = "line")
# dev.off()
#' 
#' ## Semi-parametric AFT models
#' 
#' It turns out you need to be quite careful with using time or
#' log(time) when calling different functions. Those that take
#' the linear regression model as starting point need log(time),
#' because the outcome is log(time), with the issue being that
#' this outcome sometimes is right-censored.
tic("lss from {lss2}")
ls.lss <- lss(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  mcsize = 100
)
toc()
ls.lss
#' The package aftgee uses time itself (not log(time)).
#' Computation is much faster.
tic("aftgee from {aftgee}")
aft1 <- aftgee(Surv(yrs, stat) ~ hepato + edema, data = pbc)
toc()
summary(aft1)
#' The package aftsem again needs log(time) to be specified. It
#' seems much faste than aftgee, but that is only when no SE's
#' are computed.
tic("aftsem from {aftsem}, buckley")
aftsem1 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "buckley"
)
toc()
summary(aftsem1)
#' It turns out that 15 iterations is not enough for convergence
#' (the same is true for the jin method coming up), but it is not
#' really clear how to adjust this, so I am letting this go.
tic("aftsem from {aftsem}, jin")
aftsem2 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "jin"
)
toc()
summary(aftsem2)
tic("aftsem from {aftsem}, jin, with 100 resamples")
aftsem2 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "jin", resample = 100
)
toc()
summary(aftsem2)
tic("aftsem from {aftsem}, gehan")
aftsem3 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "gehan"
)
toc()
summary(aftsem3)
tic("aftsem from {aftsem}, gehan, with 100 resamples")
aftsem3 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "gehan", resample = 100
)
toc()
summary(aftsem3)
tic("aftsem from {aftsem}, gehan-heller")
aftsem4 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "gehan-heller"
)
toc()
summary(aftsem4)
tic("aftsem from {aftsem}, gehan-poly")
aftsem5 <- aftsem(
  Surv(log(yrs), stat) ~ hepato + edema, data = pbc,
  method = "gehan-poly"
)
toc()
summary(aftsem5)
#' Finally, we have the Buckley-James estimator, implemented in
#' the rms package. This packages uses time, not log(time).
tic("bj from {rms}")
ls.bj <- bj(Surv(yrs, stat) ~ hepato + edema, data = pbc)
toc()
ls.bj

#' Here is an attempt to obtain model-based survival curves from
#' the Buckley-James estimator (obtained from ls.bj).
scaled_yrs <- pbc$yrs
S00 <- survfit(Surv(scaled_yrs, pbc$stat) ~ 1)
scaled_yrs <- pbc$yrs * exp( ls.bj$coef[3] )
S01 <- survfit(Surv(scaled_yrs, pbc$stat) ~ 1)
scaled_yrs <- pbc$yrs * exp( ls.bj$coef[2] )
S10 <- survfit(Surv(scaled_yrs, pbc$stat) ~ 1)
scaled_yrs <- pbc$yrs * exp( ls.bj$coef[2] + ls.bj$coef[3])
S11 <- survfit(Surv(scaled_yrs, pbc$stat) ~ 1)

# pdf("AFT_semipar_survival.pdf")
# Using "step" for all four curves here (the original base R code
# only used type = "s" for S00; S01/S10/S11 were drawn with
# straight lines between survfit points, inconsistent with every
# other survival curve plot in this script - likely an
# oversight).
df_bj <- combine_he(
  data.frame(time = S00$time, surv = S00$surv),
  data.frame(time = S01$time, surv = S01$surv),
  data.frame(time = S10$time, surv = S10$surv),
  data.frame(time = S11$time, surv = S11$surv)
)
p_bj <- plot_he_survival(
  df_bj, "Semi-parametric AFT", type = "step"
)
p_bj
# dev.off()

# Try out plot function of flexsurvreg
nd <- data.frame(hepato = c(0, 0, 1, 1), edema = c(0, 1, 0, 1))
# nd's row order (00, 01, 10, 11) matches he_levels, so the four
# summary() curves below can be combined directly via
# combine_he().
flexsurv_he_df <- function(fit, t = tseq) {
  s <- summary(fit, newdata = nd, t = t, ci = FALSE)
  combine_he(
    data.frame(time = s[[1]]$time, surv = s[[1]]$est),
    data.frame(time = s[[2]]$time, surv = s[[2]]$est),
    data.frame(time = s[[3]]$time, surv = s[[3]]$est),
    data.frame(time = s[[4]]$time, surv = s[[4]]$est)
  )
}
# (Set ci = TRUE in flexsurv_he_df()/summary() above, and add a
# ribbon geom in plot_he_survival(), to restore confidence
# bands.)

plot_he_survival(
  flexsurv_he_df(fitw), "Weibull AFT", type = "line"
)

plot_he_survival(
  flexsurv_he_df(fitl), "Log-logistic", type = "line"
)

plot_he_survival(
  flexsurv_he_df(fitgg), "Generalized gamma", type = "line"
)

plot_he_survival(
  flexsurv_he_df(fitn), "Log-normal", type = "line"
)

plot_he_survival(
  flexsurv_he_df(fitf), "Generalized F", type = "line"
)

plot_he_survival(
  flexsurv_he_df(fitg), "Gompertz", type = "line"
)

#' # Combined figure
#'
#' For ease of reference in the text, we also combine six of the
#' survival curve plots above - Cox, stratified Cox, the two
#' additive hazards (OLS) variants, additive hazards (MLE), and
#' semi-parametric AFT - into a single 3x2 figure, panel-
#' labelled (a)-(f) for later reference. The "Hepato / edema"
#' legend (which is identical across all six panels) is
#' collected into a single shared legend rather than repeated
#' six times.

# Prepares one panel for the combined figure: imposes the same
# x-axis breaks/range on every panel (some of the six plots
# otherwise pick different default breaks, e.g. 0/4/8/12 instead
# of 0/3/6/9/12, so 12 doesn't line up horizontally across
# panels), and adds the (a)-(f) tag directly inside the panel via
# annotate() at x = Inf, y = -Inf. This anchors the tag to the
# panel's own corner, so the tag-to-border distance is identical
# for every panel regardless of differences in title height or
# axis labels between panels - unlike patchwork's built-in
# plot_annotation(tag_levels = ...)/plot.tag.position, which
# measures position relative to the whole plot (title + axes +
# margins included) and so is not guaranteed to look the same
# distance from the panel border across differently-sized plots.
# show_ylab = FALSE additionally removes the y-axis title, for
# the right-hand column of the figure.
tag_panel <- function(p, tag, show_ylab = TRUE) {
  p <- p +
    scale_x_continuous(breaks = seq(0, 12, by = 3)) +
    coord_cartesian(xlim = c(0, 12), ylim = c(0, 1)) +
    annotate(
      "text", x = Inf, y = -Inf, label = tag,
      hjust = 1.3, vjust = -0.6, fontface = "bold"
    )
  if (!show_ylab) p <- p + labs(y = NULL)
  p
}

combined_fig <-
  (tag_panel(p_cox, "(a)") |
     tag_panel(p_scox, "(b)", show_ylab = FALSE)) /
  (tag_panel(p_ah1, "(c)") |
     tag_panel(p_ah2, "(d)", show_ylab = FALSE)) /
  (tag_panel(p_ah3, "(e)") |
     tag_panel(p_bj, "(f)", show_ylab = FALSE)) +
  plot_layout(guides = "collect")
combined_fig <- combined_fig & theme(legend.position = "bottom")
# pdf("Combined_survival_curves.pdf", width = 8, height = 10)
combined_fig
# dev.off()
