rm(list = ls())

library(shiny)
library(RColorBrewer)
library(deSolve)
library(ape)
library(phytools)
library(plotrix)
library(castor)
set.seed(123)

palette(alpha(brewer.pal(5, "Set1"), 0.75))

pulled <-
  function(birth = function(x)
    1,
    death = function(x)
      0.5,
    rho = 1,
    t = seq(0, 100, 0.1)) {
    params <- numeric(0)
    # The following function is passed to the R function ode()
    fn <- function(t, e, params) {
      with(as.list(c(e, params)),
           {
             g <- death(t) - e * (birth(t) + death(t)) + e * e * birth(t)
             return(list(g))
           })
    }
    E0 <- 1 - rho # initial value for E
    Enum <- ode(E0, t, fn, params) # Numerical integration of e
    PSRfull <- birth(t) * (1 - Enum[, 2]) # Pulled speciation rate
    # To get all vectors with the same length
    # we suppress the last value of the vector because only intervals are considered for the derivative
    PSR <- PSRfull[-length(PSRfull)]
    T <- t[-length(t)]
    dt <- T[2] - T[1] # to get the time interval
    rdlambda <- (1 / birth(T)) *  diff(birth(t)) / (dt)
    #diff(birth(t-1)) / (dt * birth(t-1)) # Relative derivative of the PSR
    PDR <-
      birth(T) - death(T) + rdlambda   # Pulled diversification rate
    PER <- birth(0) / (1 - E0) - PDR   # Pulled extinction rate
    return(list(
      "PSR" = PSR,
      "PER" = PER,
      "PDR" = PDR
    ))
  }

smooth_stepwise <- function(r0, r1, Tshift, alpha, t) {
  (r0 * exp(alpha * (t - Tshift)) + r1) / (exp(alpha * (t - Tshift)) + 1)
}

pdf("figures/fig4.pdf", height = 10, width = 8)
tmax <- 300
seq_time <- seq(0, tmax, 0.1)
spe1 <-
  function(x)
    smooth_stepwise(
      r0 = 0.04,
      r1 = 0.05,
      Tshift = 100,
      alpha = 0.5,
      x
    )
ext1 <- function(x)
  0.03 + x - x
div1 <- function(x)
  spe1(x) - ext1(x)
pulled_rates1 <- pulled(birth = spe1, death = ext1, t = seq_time)
PSR1 <- pulled_rates1$PSR
PER1 <- pulled_rates1$PER
PDR1 <- pulled_rates1$PDR

t <- -seq_time[-length(seq_time)]
par(mfrow = c(3, 2))
# True rates
Ymin <- 0.9 * min(spe1(-t), ext1(-t), div1(-t)) - 0.03
Ymax <- 1.1 * max(spe1(-t), ext1(-t), div1(-t)) + 0.02
plot(
  spe1(-t) ~ t,
  type = "l",
  col = 2,
  ylim = c(Ymin, Ymax),
  ylab = "",
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(t)), labels=rev(pretty(range(-t))))
lines(ext1(-t) ~ t, col = 1)
lines(div1(-t) ~ t, col = "black", lwd = 2)
legend(
  "topleft",
  legend = c("speciation rate",
             "extinction rate",
             "diversification rate"),
  lty = 1,
  col = c(2, 1, "black"),
  bty = 'n'
)
mtext("a)",
      side = 3,
      line = 2,
      adj = -0.1)

axis(side = 3,
     at = pretty(range(-rev((
       c(0, tmax)
     )))),
     labels = pretty(range(rev((
       c(0, tmax)
     )))))
mtext("Time (t)",
      side = 3,
      line = 2.5,
      cex = 0.65)

# True and pulled speciation rates (+ true diversification rate)
plot(
  spe1(-t) ~ t,
  type = "l",
  col = 2,
  ylim = c(Ymin, Ymax),
  ylab = "",
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(t)), labels=rev(pretty(range(-t))))
lines(PSR1 ~ t, col = 2, lty = 2)
lines(div1(-t) ~ t, col = "black", lwd = 2)
legend(
  "topleft",
  legend = c(
    "diversification rate",
    "speciation rate",
    "pulled speciation rate"
  ),
  lty = c(1, 1, 2),
  col = c("black", 2, 2),
  bty = 'n'
)
mtext("b)",
      side = 3,
      line = 2,
      adj = -0.1)


# True and pulled diversification rates
plot(
  div1(-t) ~ t,
  type = "l",
  lwd = 2,
  col = "black",
  ylim = c(Ymin, Ymax),
  ylab = "",
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(t)), labels=rev(pretty(range(-t))))
lines(PDR1 ~ t, col = "black", lty = 2)
legend(
  "topleft",
  legend = c("diversification rate", "pulled diversification rate"),
  lty = c(1, 2),
  col = "black",
  bty = 'n'
)
mtext("c)",
      side = 3,
      line = 2,
      adj = -0.1)


# True and pulled extinction rates
plot(
  ext1(-t) ~ t,
  type = "l",
  col = 1,
  ylim = c(Ymin, Ymax),
  ylab = "",
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(t)), labels=rev(pretty(range(-t))))
lines(PER1 ~ t, col = 1, lty = 2)
legend(
  "topleft",
  legend = c("extinction rate", "pulled extinction rate"),
  lty = c(1, 2),
  col = 1,
  bty = 'n'
)
mtext("d)",
      side = 3,
      line = 2,
      adj = -0.1)

# Plotss of the simulated LTT

# The function must be rewritten as a function of tmax -x to be passed to rbdtree()
spe1rev <- function(x)
  spe1(tmax - x)
ext1rev <- function(x)
  ext1(tmax - x)

# Function to count how many species are present at time t given ltt
count_species <- function(ltt, t) {
  length(which(ltt < t))
}
tmax <- 300
Nrep <- 50 # Number of simulated LTT plots
ts <- 0.01 - c(tmax:0) # x-axis to plot the mean LTT plot

#Nrep is the number of individual LTT plots, sp and deg are the span and degree paramaters for the loess function to compute the local slopes of the LTT
deg = 2
sp = 0.75

tree <-
  rbdtree(
    Tmax = tmax,
    birth = spe1rev,
    death = ext1rev,
    eps = 0.0001
  )

n <- c(1:tree$Nnode)
tsim <- -sort(branching.times(tree), decreasing = T)
meanLTT <- sapply(ts, function(x)
  count_species(tsim, x)) / Nrep

#slope of the LTT
# first loess smoothing
fit <- loess(log(n) ~ tsim, span = sp, degree = deg)
fitted_logn <- predict(fit, tsim)
slope <- diff(fitted_logn) / diff(tsim)
l <- length(slope)
tslope <- tsim[-1] # to keep the same length as slope
numsim <- rep(1, l)


plot(
  c(n, max(n)) ~ c(tsim, 0),
  log = "y",
  xlim = c(min(ts), max(ts)),
  ylim = c(1, 10000),
  type = "l",
  ylab = "log # of species",
  col = "grey",
  cex = 0.5,
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
)
axis(1, at=pretty(range(t)), labels=rev(pretty(range(-t))))


for (i in 1:(Nrep - 1)) {
  meanLTT <-
    meanLTT + sapply(ts, function(x)
      count_species(tsim, x)) / Nrep
  tree <-
    rbdtree(
      Tmax = tmax,
      birth = spe1rev,
      death = ext1rev,
      eps = 0.001
    )
  tsim <- -sort(branching.times(tree), decreasing = T)
  n <- c(1:tree$Nnode)
  lines(c(n, max(n)) ~ c(tsim, 0), col = "grey")
  fit <- loess(log(n) ~ tsim, span = sp, degree = deg)
  fitted_logn <- predict(fit, tsim)
  new_slope <- diff(fitted_logn) / diff(tsim)
  l <- length(new_slope)
  slope <- c(slope, new_slope)
  tslope <- c(tslope, tsim[-1])
  numsim <- c(numsim, rep(i + 1, l))
}

lines(meanLTT ~ ts , lwd = 2)
legend(
  "topleft",
  legend = c(
    "simulated LTT",
    "mean LTT"
  ),
  lty = c(1, 1),
  col = c('grey','black'),
  bty = 'n'
)

mtext("e)",
      side = 3,
      line = 2,
      adj = -0.1)
# Plot of the slope of the LTT

sim_slope <- data.frame(cbind(numsim, tslope, slope))

plot(
  NULL,
  NULL,
  xlim = c(-tmax, 0),
  ylab = "Slope of LTT",
  ylim = c(Ymin, Ymax),
  xlab = bquote("Age (" * tau * ")"),
  xaxt = "n"
) # Pulled-speciation rate
for (i in c(1:Nrep))
  points(sim_slope[sim_slope$numsim == i, ]$slope ~ sim_slope[sim_slope$numsim ==
                                                                i, ]$tslope, col = alpha("grey",0.5),pch=16)
axis(1, at = pretty(range(t)), labels = rev(pretty(range(-t))))
lines(-seq_time[-1],
      PSR1,
      col = 2,
      lty = 2,
      lwd = 1) # Pulled-speciation rate
lines(ts, div1(-ts), lty = 1, lwd = 1) # Diversification rate
lines(ts,
      spe1(-ts),
      col = 2,
      lty = 1,
      lwd = 1) # Speciation rate
legend(
  "topleft",
  legend = c(
    "diversification rate",
    "speciation rate",
    "pulled speciation rate",
    "LTT slope"
  ),
  lty = c(1, 1, 2, NA),
  col = c("black", 2, 2,alpha('grey',0.5)),
  pch= c(NA,NA,NA,16),
  bty = 'n'
)
mtext("f)",
      side = 3,
      line = 2,
      adj = -0.1)

dev.off()
