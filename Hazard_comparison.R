tseq <- seq(0, 5, by = 0.01)

h0 <- function(t) 0.1 * (t - 1)^2 + 0.1
H0 <- function(t) 0.1/3 * (t - 1)^3 + 0.1 * t + 0.1/3
S0 <- function(t) exp(-H0(t))

plot(tseq, h0(tseq), type = "l", lwd = 2)
plot(tseq, H0(tseq), type = "l", lwd = 2)
plot(tseq, S0(tseq), type = "l", lwd = 2)

RMST5 <- integrate(S0, lower = 0, upper = 5)
RMST5
RMST5 <- RMST5$value
## Find beta such that RMST5 for X=1 decreases with 1
# Cox model
f <- function(bet) {
  S1_cox <- function(t) exp(-H0(t) * exp(bet))
  integrate(S1_cox, lower = 0, upper = 5)$value - RMST5 + 1
}
ur <- uniroot(f, lower = 0, upper = 1)
ur
bet_cox <- ur$root
hplus_cox <- function(t) h0(t) * exp(bet_cox)
hminus_cox <- function(t) h0(t) * exp(-bet_cox)
Splus_cox <- function(t) S0(t) ^ exp(bet_cox)
Sminus_cox <- function(t) S0(t) ^ exp(-bet_cox)
Hplus_cox <- function(t) H0(t) * exp(bet_cox)
Hminus_cox <- function(t) H0(t) * exp(-bet_cox)
# Additive hazards model
f <- function(bet) {
  S1_ah <- function(t) exp(-H0(t) - bet * t)
  integrate(S1_ah, lower = 0, upper = 5)$value - RMST5 + 1
}
ur <- uniroot(f, lower = 0, upper = 1)
ur
bet_ah <- ur$root
hplus_ah <- function(t) h0(t) + bet_ah
hminus_ah <- function(t) h0(t) - bet_ah
Splus_ah <- function(t) S0(t) * exp(-bet_ah * t)
Sminus_ah <- function(t) S0(t) * exp(bet_ah * t)
Hplus_ah <- function(t) H0(t) + bet_ah * t
Hminus_ah <- function(t) H0(t) - bet_ah * t
# AFT model
f <- function(bet) {
  S1_AFT <- function(t) S0(t * exp(-bet))
  integrate(S1_AFT, lower = 0, upper = 5)$value - RMST5 + 1
}
ur <- uniroot(f, lower = -1, upper = 0)
ur
bet_AFT <- ur$root
hplus_AFT <- function(t) h0(t * exp(-bet_AFT)) * exp(-bet_AFT)
hminus_AFT <- function(t) h0(t * exp(bet_AFT)) * exp(bet_AFT)
Splus_AFT <- function(t) S0(t * exp(-bet_AFT))
Sminus_AFT <- function(t) S0(t * exp(bet_AFT))
Hplus_AFT <- function(t) H0(t * exp(-bet_AFT))
Hminus_AFT <- function(t) H0(t * exp(bet_AFT))

### Plot
tseq <- seq(0, 5, by = 0.01)

# pdf("Survival_curves.pdf", width = 8, height = 3)
## Survival curves
my_ylim = c(0, 1)
par(mfrow = c(1, 3))
# Cox
plot(tseq, S0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Splus_cox(tseq), lty = 3)
lines(tseq, Sminus_cox(tseq), lty = 2)
title(main = "Proportional hazards")
# Additive hazards
plot(tseq, S0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Splus_ah(tseq), lty = 3)
lines(tseq, Sminus_ah(tseq), lty = 2)
title(main = "Additive hazards")
# AFT
plot(tseq, S0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Splus_AFT(tseq), lty = 3)
lines(tseq, Sminus_AFT(tseq), lty = 2)
title(main = "Accelerated failure time")
# dev.off()

# pdf("Hazard_curves.pdf", width = 8, height = 3)
## Hazards (original scale)
my_ylim = c(-0.2, 5)
par(mfrow = c(1, 3))
# Cox
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, hplus_cox(tseq), lty = 3)
lines(tseq, hminus_cox(tseq), lty = 2)
title(main = "Proportional hazards")
# Additive hazards
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, hplus_ah(tseq), lty = 3)
lines(tseq, hminus_ah(tseq), lty = 2)
title(main = "Additive hazards")
# AFT
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, hplus_AFT(tseq), lty = 3)
lines(tseq, hminus_AFT(tseq), lty = 2)
title(main = "Accelerated failure time")
# dev.off()

# pdf("Hazard_curves_log.pdf", width = 8, height = 3)
## Hazards (log-scale y-axis)
my_ylim = c(0.025, 5)
par(mfrow = c(1, 3))
# Cox
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function", log = "y")
lines(tseq, hplus_cox(tseq), lty = 3)
lines(tseq, hminus_cox(tseq), lty = 2)
title(main = "Proportional hazards")
# Additive hazards
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function", log = "y")
lines(tseq, hplus_ah(tseq), lty = 3)
lines(tseq, hminus_ah(tseq), lty = 2)
title(main = "Additive hazards")
# AFT
plot(tseq, h0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function", log = "y")
lines(tseq, hplus_AFT(tseq), lty = 3)
lines(tseq, hminus_AFT(tseq), lty = 2)
title(main = "Accelerated failure time")
# dev.off()

# pdf("Cum_hazard_curves.pdf", width = 8, height = 3)
## Hazards (original scale)
my_ylim = c(-0.2, 5)
par(mfrow = c(1, 3))
# Cox
plot(tseq, H0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Hplus_cox(tseq), lty = 3)
lines(tseq, Hminus_cox(tseq), lty = 2)
title(main = "Proportional hazards")
# Additive hazards
plot(tseq, H0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Hplus_ah(tseq), lty = 3)
lines(tseq, Hminus_ah(tseq), lty = 2)
title(main = "Additive hazards")
# AFT
plot(tseq, H0(tseq), type = "l", lwd = 2, ylim = my_ylim,
     xlab = "Time", ylab = "Hazard function")
lines(tseq, Hplus_AFT(tseq), lty = 3)
lines(tseq, Hminus_AFT(tseq), lty = 2)
title(main = "Accelerated failure time")
# dev.off()

# When is the hazard for additive hazards model for X=-1 negative?
hminus_ah(tseq[223:224])
tseq[223:224]
# When is the survival curve for additive hazards model for X=-1 below one?
Sminus_ah(tseq[343:344])
tseq[343:344]

# Wat ik nog extra kan doen (misschien alleen voor github?) is het volgende.
# Ik kan de meer flexibele \ah en \aft modellen gebruiken, en ervoor zorgen
# dat ik in beide gevallen (\eqref{eq:ah} voor \ah, \eqref{eq:td_AFT} voor AFT)
# $\beta(t)$ zodanig kies dat zowel voor $X=0$ als $X=1$ de hazards en dus
# survival curves exact gelijk zijn voor de drie modellen. Ik hoop dat ik dat
# voor elkaar krijg. Voor \ah kan dat makkelijk, voor AFT weet ik nog niet.
# Vervolgens laten zien hoe de (log) hazard en survival curves eruit zien
# voor $X=-1$. Heb ik nog niet gedaan, maar kan ook inzicht opleveren.

