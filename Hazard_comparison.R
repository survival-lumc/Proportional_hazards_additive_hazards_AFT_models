#####
# Comparison of hazard, cumulative hazard, and survival functions
# across three commonly used survival models:
#   - Proportional hazards (Cox)
#   - Additive hazards
#   - Accelerated failure time (AFT)
#
# A common baseline hazard h0(t) is used for all models, and for
# each model the covariate effect (X = +1 / X = -1, relative to
# baseline) is calibrated so that the restricted mean survival
# time (RMST) up to tau = 5 changes by exactly 1 unit for X = +1.
# This lets the three models be compared on a common footing.
#
# The script then illustrates non-collapsibility: a discrete
# frailty variable U is introduced, and population-averaged
# ("marginal") hazards are derived by integrating out U from the
# individual-level ("conditional") hazards. The script shows how,
# in contrast to the additive hazards and AFT models, the
# population-averaged hazard ratio for the Cox model changes over
# time even when the individual-level (conditional) hazard ratio
# is constant - the well-known non-collapsibility property of the
# hazard ratio.
#####

library(tidyverse)
library(patchwork)
library(numDeriv)

#####
# Baseline functions
#####

tseq <- seq(0, 5, by = 0.001)

h0 <- function(t) 0.1 * (t - 1)^2 + 0.1
H0 <- function(t) 0.1/3 * (t - 1)^3 + 0.1 * t + 0.1/3
S0 <- function(t) exp(-H0(t))

plot(tseq, h0(tseq), type = "l", lwd = 2, main = "Hazard")
plot(
  tseq, H0(tseq), type = "l", lwd = 2,
  main = "Cumulative hazard"
)
plot(tseq, S0(tseq), type = "l", lwd = 2, main = "Survival",
     ylim = c(0, 1))

# Baseline RMST until tau = 5
RMST5 <- integrate(S0, lower = 0, upper = 5)
RMST5
RMST5 <- RMST5$value
## Find beta such that RMST5 for X=1 decreases with 1
# For each model, beta is found by root-finding (uniroot) so that
# the RMST up to tau = 5 for X = +1 is exactly 1 unit lower than
# the baseline RMST. This calibrates the three models to a
# comparable effect size. The search interval differs by model
# because of how beta enters the survival function: for Cox and
# additive hazards, increasing beta from 0 upward decreases
# survival (so the root lies in [0, 1]), while for AFT, survival
# decreases as beta decreases below 0 (time is "stretched" by
# exp(-beta), so a negative beta speeds up the time scale and
# lowers survival) - hence the root lies in [-1, 0].
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
# AFT model (Accelerated Failure Time): covariate X rescales time
# itself rather than the hazard, so S1(t) = S0(t * exp(-bet)) for
# the AFT model.
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
# Prepare data
plot_data <- tibble(t = tseq) |>
  mutate(
    S0 = S0(t),
    # Cox / Proportional hazards
    Splus_cox = Splus_cox(t),
    Sminus_cox = Sminus_cox(t),
    # Additive hazards
    Splus_ah = Splus_ah(t),
    Sminus_ah = Sminus_ah(t),
    # AFT
    Splus_AFT = Splus_AFT(t),
    Sminus_AFT = Sminus_AFT(t)
  )

# Common theme
theme_surv <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Proportional hazards (Cox)
p_cox <- plot_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = S0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Splus_cox, linetype = "X = +1")) +
  geom_line(aes(y = Sminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Survival function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# Additive hazards
p_ah <- plot_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = S0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Splus_ah, linetype = "X = +1")) +
  geom_line(aes(y = Sminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# AFT
p_aft <- plot_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = S0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Splus_AFT, linetype = "X = +1")) +
  geom_line(aes(y = Sminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# Combine with patchwork
p_cox + p_ah + p_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save as PDF (optional)
# ggsave("Survival_curves.pdf", width = 10, height = 4)

### Hazard curves
# Prepare data for hazard curves
hazard_data <- tibble(t = tseq) |>
  mutate(
    h0 = h0(t),
    # Cox / Proportional hazards
    hplus_cox = hplus_cox(t),
    hminus_cox = hminus_cox(t),
    # Additive hazards
    hplus_ah = hplus_ah(t),
    hminus_ah = hminus_ah(t),
    # AFT
    hplus_AFT = hplus_AFT(t),
    hminus_AFT = hminus_AFT(t)
  )

# Proportional hazards (Cox)
p_haz_cox <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_cox, linetype = "X = +1")) +
  geom_line(aes(y = hminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Additive hazards
p_haz_ah <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_ah, linetype = "X = +1")) +
  geom_line(aes(y = hminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# AFT
p_haz_aft <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_AFT, linetype = "X = +1")) +
  geom_line(aes(y = hminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Combine panels
p_haz_cox + p_haz_ah + p_haz_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save as PDF (optional)
# ggsave("Hazard_curves.pdf", width = 10, height = 4)

### Hazard curves on the log scale
# Prepare data for hazard curves (reuse hazard_data from before)

# Proportional hazards (Cox) - log scale
p_haz_log_cox <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_cox, linetype = "X = +1")) +
  geom_line(aes(y = hminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = c(0.025, 5)) +
  labs(
    x = "Time",
    y = "Hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  theme_surv

# Additive hazards - log scale
p_haz_log_ah <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_ah, linetype = "X = +1")) +
  geom_line(aes(y = hminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = c(0.025, 5)) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  theme_surv

# AFT - log scale
p_haz_log_aft <- hazard_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = h0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = hplus_AFT, linetype = "X = +1")) +
  geom_line(aes(y = hminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = c(0.025, 5)) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  theme_surv

# Combine panels
p_haz_log_cox + p_haz_log_ah + p_haz_log_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save as PDF (optional)
# ggsave("Hazard_curves_log.pdf", width = 10, height = 4)

### Cumulative hazards plots
# Prepare data for cumulative hazard curves
cumhaz_data <- tibble(t = tseq) |>
  mutate(
    H0 = H0(t),
    # Cox / Proportional hazards
    Hplus_cox = Hplus_cox(t),
    Hminus_cox = Hminus_cox(t),
    # Additive hazards
    Hplus_ah = Hplus_ah(t),
    Hminus_ah = Hminus_ah(t),
    # AFT
    Hplus_AFT = Hplus_AFT(t),
    Hminus_AFT = Hminus_AFT(t)
  )

# Proportional hazards (Cox)
p_cumhaz_cox <- cumhaz_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = H0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Hplus_cox, linetype = "X = +1")) +
  geom_line(aes(y = Hminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Cumulative hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Additive hazards
p_cumhaz_ah <- cumhaz_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = H0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Hplus_ah, linetype = "X = +1")) +
  geom_line(aes(y = Hminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# AFT
p_cumhaz_aft <- cumhaz_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = H0, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Hplus_AFT, linetype = "X = +1")) +
  geom_line(aes(y = Hminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Combine panels
p_cumhaz_cox + p_cumhaz_ah + p_cumhaz_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save as PDF (optional)
# ggsave("Cum_hazard_curves.pdf", width = 10, height = 4)

# When is the hazard for additive hazards model for X=-1
# negative? Indices found by manual inspection of hminus_ah(tseq)
# to bracket the point where it crosses zero (tseq <- seq(0, 5,
# by = 0.001), so index i corresponds to time tseq[i] = (i - 1) *
# 0.001, i.e. t ~ 2.22-2.23 here).
hminus_ah(tseq[c(2221, 2231)])
tseq[c(2221, 2231)]
# When is the survival curve for additive hazards model for X=-1
# below one? Same approach: indices bracket the relevant crossing
# point (t ~ 3.42-3.43), found by manual inspection of
# Sminus_ah(tseq).
Sminus_ah(tseq[c(3421, 3431)])
tseq[c(3421, 3431)]

######
# Illustration of non-collapsibility
######
# A discrete frailty variable U is mixed into the population to
# create heterogeneity that is *not* observed/conditioned on.
# Averaging ("marginalizing") the individual-level
# (U-conditional) hazards over U gives population-averaged
# ("marginal") hazards. For the Cox model in particular, the
# population-averaged hazard ratio for X is *not* constant over
# time even though the individual-level (U-conditional) hazard
# ratio is - this is the non-collapsibility of the Cox hazard
# ratio. This contrasts with the additive hazards and AFT models,
# where the corresponding population-averaged effect measures
# stay constant.

# Frailty variable U takes values -1, 0, 1 with probabilities
# pmin1, p0, p1
pmin1 <- p1 <- 0.5
p0 <- 0

# Take effect of U the same as that of X (could be different)
gam_cox <- bet_cox
gam_ah <- bet_ah
gam_AFT <- bet_AFT

# Hazard ratios (r for ratio) related to values of U
rmin1 <- exp(-gam_cox)
r0 <- 1
r1 <- exp(gam_cox)

# And later useful, pr for product of p and r
prmin1 <- pmin1 * rmin1
pr0 <- p0 * r0
pr1 <- p1 * r1

# Laplace transform of the frailty distribution U, evaluated at
# the cumulative baseline hazard, gives the population-averaged
# survival function under the Cox model: S_av(t) =
# E_U[exp(-U-weighted H(t))] = sum_u P(U=u) * exp(-r_u * H(t)) =
# Laplace(H(t)). Laplace_der is its derivative with respect to
# its argument c, and Laplace_ratio = -Laplace_der(c) /
# Laplace(c) is the multiplicative factor that converts the
# conditional (individual-level) hazard into the
# population-averaged (marginal) hazard via hav(t) = h(t) *
# Laplace_ratio(H(t)).
Laplace <- function(c) pmin1 * exp(-rmin1 * c) +
  p0 * exp(-r0 * c) + p1 * exp(-r1 * c)
Laplace_der <- function(c) -(prmin1 * exp(-rmin1 * c) +
                               pr0 * exp(-r0 * c) + pr1 * exp(-r1 * c))
Laplace_ratio <- function(c) -Laplace_der(c) / Laplace(c)

# Cox model
Sav_cox <- function(t) Laplace(H0(t))
Savplus_cox <- function(t) Laplace(Hplus_cox(t))
Savminus_cox <- function(t) Laplace(Hminus_cox(t))
Hav_cox <- function(t) -log(Sav_cox(t))
Havplus_cox <- function(t) -log(Savplus_cox(t))
Havminus_cox <- function(t) -log(Savminus_cox(t))
hav_cox <- function(t) h0(t) * Laplace_ratio(H0(t))
havplus_cox <- function(t) {
  hplus_cox(t) * Laplace_ratio(Hplus_cox(t))
}
havminus_cox <- function(t) {
  hminus_cox(t) * Laplace_ratio(Hminus_cox(t))
}

# Additive hazards model Renamed to *_U to avoid overwriting
# Splus_ah/Sminus_ah defined earlier (those represent the effect
# of X; these represent the effect of U). Currently identical in
# value since gam_ah <- bet_ah above, but kept as separate
# objects so the two effects can be set independently.
Splus_ah_U <- function(t) S0(t) * exp(-gam_ah * t)
Sminus_ah_U <- function(t) S0(t) * exp(gam_ah * t)
Sav_ah <- function(t) {
  pmin1 * Sminus_ah_U(t) + p0 * S0(t) + p1 * Splus_ah_U(t)
}
Savplus_ah <- function(t) Sav_ah(t) * exp(-bet_ah * t)
Savminus_ah <- function(t) Sav_ah(t) * exp(bet_ah * t)
Hav_ah <- function(t) -log(Sav_ah(t))
Havplus_ah <- function(t) Hav_ah(t) + bet_ah * t
Havminus_ah <- function(t) Hav_ah(t) - bet_ah * t
# hav_af <- sapply(tseq, function(t) numDeriv::grad(Hav_ah, t))
# plot(tseq, h0(tseq), type = "l")
# lines(tseq, hav_af, type = "l", col = 2)
hav_ah <- function(t) {
  den <- pmin1 * exp(gam_ah * t) + p0 + p1 * exp(-gam_ah * t)
  num <- pmin1 * gam_ah * exp(gam_ah * t) -
    p1 * gam_ah * exp(-gam_ah * t)
  return(h0(t) - num / den)
}
havplus_ah <- function(t) hav_ah(t) + bet_ah
havminus_ah <- function(t) hav_ah(t) - bet_ah
# plot(tseq, h0(tseq), type = "l")
# lines(tseq, hav_af(tseq), type = "l", col = 2)

# AFT model
Sp_AFT <- function(t) S0(t * exp(-gam_AFT))
Sm_AFT <- function(t) S0(t * exp(gam_AFT))
Sav_AFT <- function(t) {
  pmin1 * Sm_AFT(t) + p0 * S0(t) + p1 * Sp_AFT(t)
}
Spplus_AFT <- function(t) S0(t * exp(-bet_AFT -gam_AFT))
S0plus_AFT <- function(t) S0(t * exp(-bet_AFT))
Smplus_AFT <- function(t) S0(t * exp(-bet_AFT + gam_AFT))
Savplus_AFT <- function(t) {
  pmin1 * Smplus_AFT(t) + p0 * S0plus_AFT(t) +
    p1 * Spplus_AFT(t)
}
Spminus_AFT <- function(t) S0(t * exp(bet_AFT -gam_AFT))
S0minus_AFT <- function(t) S0(t * exp(bet_AFT))
Smminus_AFT <- function(t) S0(t * exp(bet_AFT + gam_AFT))
Savminus_AFT <- function(t) {
  pmin1 * Smminus_AFT(t) + p0 * S0minus_AFT(t) +
    p1 * Spminus_AFT(t)
}
Hav_AFT <- function(t) -log(Sav_AFT(t))
Havplus_AFT <- function(t) -log(Savplus_AFT(t))
Havminus_AFT <- function(t) -log(Savminus_AFT(t))
hav_AFT <- sapply(
  tseq, function(t) numDeriv::grad(Hav_AFT, t)
)
havplus_AFT <- sapply(
  tseq, function(t) numDeriv::grad(Havplus_AFT, t)
)
havminus_AFT <- sapply(
  tseq, function(t) numDeriv::grad(Havminus_AFT, t)
)
# plot(tseq, h0(tseq), type = "l")
# lines(tseq, hav_AFT, type = "l", col = 2)

# Prepare data for averaged survival curves
surv_av_data <- tibble(t = tseq) |>
  mutate(
    # Cox / Proportional hazards
    Sav_cox = Sav_cox(t),
    Savplus_cox = Savplus_cox(t),
    Savminus_cox = Savminus_cox(t),
    # Additive hazards
    Sav_ah = Sav_ah(t),
    Savplus_ah = Savplus_ah(t),
    Savminus_ah = Savminus_ah(t),
    # AFT
    Sav_AFT = Sav_AFT(t),
    Savplus_AFT = Savplus_AFT(t),
    Savminus_AFT = Savminus_AFT(t)
  )

# Proportional hazards (Cox)
p_surv_av_cox <- surv_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = Sav_cox, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Savplus_cox, linetype = "X = +1")) +
  geom_line(aes(y = Savminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Survival function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# Additive hazards
p_surv_av_ah <- surv_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = Sav_ah, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Savplus_ah, linetype = "X = +1")) +
  geom_line(aes(y = Savminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# AFT
p_surv_av_aft <- surv_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = Sav_AFT, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Savplus_AFT, linetype = "X = +1")) +
  geom_line(aes(y = Savminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv

# Combine panels
p_surv_av_cox + p_surv_av_ah + p_surv_av_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save as PDF (optional)
# ggsave("Survival_curves_av.pdf", width = 10, height = 4)

### Averaged hazard curves
# Prepare data for averaged hazard curves Note: hav_AFT,
# havplus_AFT, havminus_AFT were precomputed above as plain
# numeric vectors over tseq (via numDeriv::grad), not as
# functions like their Cox/additive-hazards counterparts - they
# can only be evaluated at the tseq grid used here, not at
# arbitrary new time points.
haz_av_data <- tibble(t = tseq) |>
  mutate(
    # Cox / Proportional hazards
    hav_cox     = hav_cox(t),
    havplus_cox  = havplus_cox(t),
    havminus_cox = havminus_cox(t),
    # Additive hazards
    hav_ah      = hav_ah(t),
    havplus_ah   = havplus_ah(t),
    havminus_ah  = havminus_ah(t),
    # AFT (vectors, not functions)
    hav_AFT     = hav_AFT,
    havplus_AFT  = havplus_AFT,
    havminus_AFT = havminus_AFT
  )

# Proportional hazards (Cox) Note: y-axis range here (-0.2, 2) is
# narrower than the individual-level hazard plots above (-0.2,
# 5), since population-averaged hazards have a smaller range than
# the individual-level (conditional) hazards.
p_haz_av_cox <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = hav_cox, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_cox,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 2) +
  theme_surv

# Additive hazards
p_haz_av_ah <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = hav_ah, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_ah,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 2) +
  theme_surv

# AFT
p_haz_av_aft <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = hav_AFT, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_AFT,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(-0.2, 2) +
  theme_surv

# Combine panels
p_haz_av_cox + p_haz_av_ah + p_haz_av_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ggsave("Hazard_curves_av.pdf", width = 10, height = 4)

# Plot only of Cox model partly conditional hazards
# Note: "partly conditional" here refers to the same
# population-averaged/marginal hazard (hav_cox etc.) computed in
# the "Averaged hazard curves" section above - just re-plotted as
# a single panel with a different legend style. Renamed to
# p_haz_partly_cond_cox (rather than reusing p_haz_cox) so it
# doesn't overwrite the individual-level Cox hazard plot defined
# earlier in the script.
haz_cox_data <- tibble(t = tseq) |>
  mutate(
    hav_cox      = hav_cox(t),
    havplus_cox  = havplus_cox(t),
    havminus_cox = havminus_cox(t)
  )

p_haz_partly_cond_cox <- haz_cox_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = hav_cox, linetype = "X = 0"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_cox,  linetype = "X = 1")) +
  geom_line(aes(y = havminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "X = 0" = "solid",
      "X = 1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Hazard function",
    title = "Partly conditional hazards\nCox model",
    linetype = NULL
  ) +
  ylim(0, 1) +
  theme_surv +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1)
  )

p_haz_partly_cond_cox

# ggsave("Partly_cond_haz_Cox.pdf", width = 6, height = 6)

# Plot of hazard ratios
# Hazard ratio exp(bet_cox):
exp(bet_cox)

hr_cox_data <- tibble(t = tseq) |>
  mutate(
    hr_plus  = havplus_cox(t) / hav_cox(t),
    hr_minus = hav_cox(t) / havminus_cox(t)
  )

p_hr_cox <- hr_cox_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = hr_plus, colour = "X = 1 compared to X = 0"),
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = hr_minus, colour = "X = 0 compared to X = -1"),
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "X = 1 compared to X = 0" = "black",
      "X = 0 compared to X = -1" = "red"
    )
  ) +
  labs(
    x = "Time",
    y = "Hazard ratio",
    title = "Hazard ratios",
    colour = NULL
  ) +
  theme_surv +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1)
  )

p_hr_cox

# ggsave("Partly_cond_HR.pdf", width = 6, height = 6)

# Check for additive hazards
# looks like same number throughout
head(havplus_ah(tseq) - hav_ah(tseq))
bet_ah # namely the conditional effect of X
all.equal(havplus_ah(tseq), hav_ah(tseq) + bet_ah) # check
all.equal(havminus_ah(tseq), hav_ah(tseq) - bet_ah) # check

# It is very hard to judge the AFT nature from the hazards
# Easiest to show that these average survival curves for AFT
# still follows AFT with bet_AFT, by looking up the median (or
# any other quantile)
tmed <- min(tseq[Sav_AFT(tseq) < 0.5])
tmed
Sav_AFT(tmed)
Savplus_AFT(tmed * exp(bet_AFT))
Savminus_AFT(tmed * exp(-bet_AFT))

my_ylim_log <- c(0.025, 5)

# Renamed to p_haz_log_av_* (rather than reusing
# p_haz_log_cox/ah/aft) so these don't overwrite the
# individual-level log-hazard plots defined earlier in the
# "Hazard curves on the log scale" section.

# Proportional hazards (Cox)
p_haz_log_av_cox <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = hav_cox, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_cox,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = my_ylim_log) +
  labs(
    x = "Time",
    y = "Hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  theme_surv

# Additive hazards
p_haz_log_av_ah <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = hav_ah, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_ah,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = my_ylim_log) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  theme_surv

# AFT
p_haz_log_av_aft <- haz_av_data |>
  ggplot(aes(x = t)) +
  geom_line(
    aes(y = hav_AFT, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = havplus_AFT,  linetype = "X = +1")) +
  geom_line(aes(y = havminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  scale_y_log10(limits = my_ylim_log) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  theme_surv

# Combine panels
p_haz_log_av_cox + p_haz_log_av_ah + p_haz_log_av_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ggsave("Hazard_curves_log_av.pdf", width = 10, height = 4)

### Averaged cumulative hazard curves
# Prepare data for averaged cumulative hazard curves
cum_haz_av_data <- tibble(t = tseq) |>
  mutate(
    # Cox / Proportional hazards
    Hav_cox      = Hav_cox(t),
    Havplus_cox  = Havplus_cox(t),
    Havminus_cox = Havminus_cox(t),
    # Additive hazards
    Hav_ah       = Hav_ah(t),
    Havplus_ah   = Havplus_ah(t),
    Havminus_ah  = Havminus_ah(t),
    # AFT
    Hav_AFT      = Hav_AFT(t),
    Havplus_AFT  = Havplus_AFT(t),
    Havminus_AFT = Havminus_AFT(t)
  )

# Proportional hazards (Cox)
p_cum_haz_cox <- cum_haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = Hav_cox, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Havplus_cox,  linetype = "X = +1")) +
  geom_line(aes(y = Havminus_cox, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = "Cumulative hazard function",
    title = "Proportional hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Additive hazards
p_cum_haz_ah <- cum_haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = Hav_ah, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Havplus_ah,  linetype = "X = +1")) +
  geom_line(aes(y = Havminus_ah, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Additive hazards",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# AFT
p_cum_haz_aft <- cum_haz_av_data |>
  ggplot(aes(x = t)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_line(
    aes(y = Hav_AFT, linetype = "Baseline"),
    linewidth = 0.8
  ) +
  geom_line(aes(y = Havplus_AFT,  linetype = "X = +1")) +
  geom_line(aes(y = Havminus_AFT, linetype = "X = -1")) +
  scale_linetype_manual(
    values = c(
      "Baseline" = "solid",
      "X = +1" = "dotted",
      "X = -1" = "dashed"
    )
  ) +
  labs(
    x = "Time",
    y = NULL,
    title = "Accelerated failure time",
    linetype = NULL
  ) +
  ylim(-0.2, 5) +
  theme_surv

# Combine panels
p_cum_haz_cox + p_cum_haz_ah + p_cum_haz_aft +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ggsave("Cum_hazard_curves_av.pdf", width = 10, height = 4)