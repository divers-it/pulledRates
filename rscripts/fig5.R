rm(list = ls())
library(deSolve)
library(RColorBrewer)
palette(brewer.pal(n = 5, name = "Set1"))

# Birth and death function as a function of time: the functions you can play with

#starting value + ...

#exp( width of peak * ... ) smaller = wider

#to store birth and death functions
bl <- list()
dl <- list()

#x-20 = time of div rate change
# ^2 shape of curve

###
#recent radiation
###

birth <- function(x)
  0.3 + exp(-0.1 * (x) ^ 2)

#constant extinction
death <- function(x)
  0.05 + (x - x)

bl[[1]] <- birth
dl[[1]] <- death

###
# mass extinction
###

birth <- function(x)
  0.3 + (x - x)

death <- function(x)
  0.05 + exp(-5 * (x - 40) ^ 2)

bl[[2]] <- birth
dl[[2]] <- death

###
# Gradual increase
###

birth <- function(x)
  0.3 + exp(-0.0015 * (x) ^ 2)

death <- function(x)
  0.05 + exp(-0.0015 * (x) ^ 2)

bl[[3]] <- birth
dl[[3]] <- death


###
# Crazy
###

birth <- function(x)
  0.3 * (2.5 + exp(-0.05 * (x - 20) ^ 4) + sin(x / 7) ^ 3)
death <- function(x)
  0.28 * (2.5 + exp(-0.05 * (x - 20) ^ 4) + sin(x / 7) ^ 3)
#death <- function(x)
#  0.05 * (1.5 + 0.05*x + cos(x/20)^4 * cos(x/3)+ cos(x/5)^3)


bl[[4]] <- birth
dl[[4]] <- death

xmax = c(20, 60, 100, 100)
xmin = c(0, 20, 0, 0)
lab <- c("a)", "b)", "c)", "d)")
scenario <- c(
  "recent radiation",
  "mass extinction",
  "gradual turnover increase",
  "rapid fluctuations"
)


#####
# plot
#####

pdf("figures/fig5.pdf",height=9,width=11)

par(mfrow = c(2, 2))

for (i in 1:4) {
  birth <- bl[[i]]
  death <- dl[[i]]

  ## Solving the equation for E(t) numerically
  params <- numeric(0)
  rho = 1
  dt <- 0.1
  t <- seq(0, 100, dt) # sequence of time for integration

  # The function to give to the ODE solver
  fn <- function(t, e, params) {
    with(as.list(c(e, params)),
         {
           g <- death(t) - e * (birth(t) + death(t)) + e * e * birth(t)
           return(list(g))
         })
  }




  # initial value for e0 = 1 - rho (sampling fraction)
  E0 <- 1 - rho # initial value for E
  Enum <- ode(E0, t, fn, params) # Numerical integration of e
  PSRfull <- birth(t) * (1 - Enum[, 2]) # Pulled speciation rate

  PSR <- PSRfull[-length(PSRfull)]
  T <- t[-length(t)]
  dt <- T[2] - T[1] # to get the time interval
  rdlambda <- ( 1 / birth(T) ) *  diff(birth(t)) / ( dt )
  #diff(birth(t-1)) / (dt * birth(t-1)) # Relative derivative of the PSR
  PDR <- birth(T) - death(T) + rdlambda   # Pulled diversification rate

  divers <- birth(t) - death(t)

  PER <- rep(birth(0),length(rdlambda)) - divers[-length(divers)] - rdlambda # Pulled extinction rate

  PER <- birth(0) / (1 - E0) - PDR
  #plot

  plot(
    rev(PDR) ~ rev(-T),
    type = "l",
    lty = 2,
    ylim = c(-1, 1.5),
    xlim = -rev(c(xmin[i],xmax[i])),
    ylab = 'Rate',
    xlab = bquote("Age (" * tau * ")"),
    xaxt = "n",
    col="black"
  )
  axis(1, at=pretty(range(-rev(c(xmin[i],xmax[i])))), labels=rev(pretty(range(rev(c(xmin[i],xmax[i]))))))
  #div rate
  lines(rev(divers) ~ rev(-t),col="black")
  lines(rev(birth(T)) ~ rev(-T), col = 2)
  lines(rev(death(T)) ~ rev(-T), col = 1)
  lines(rev(PSR) ~ rev(-T), col = 2, lty = 2)
  lines(rev(PER) ~ rev(-T), col = 1, lty = 2)

  if (i == 2) {
    legend(
      "topleft",
      legend = c(
        "diversification rate",
        "pulled diversification rate",
        "speciation rate",
        "pulled speciation rate",
        "extinction rate",
        "pulled extinction rate"
      ),
      col = c("black", "black", 2, 2, 1, 1),
      lty = c(1, 2),
      bty = "n"
    )
  }

  mtext(lab[i],
        side = 3,
        line = 1,
        adj = -0.1)

  if (i == 1) {
    text(x = -2, y = 1.45, labels = scenario[i])
    axis(side=3, at = pretty(range(-rev((c(xmin[i],xmax[i]))))),labels=pretty(range(rev((c(xmin[i],xmax[i]))))))
    mtext('Time (t)', side=3, line=2.5,cex=0.8)
  }

  if (i == 2) {
    text(x = -24, y = 1.45, labels = scenario[i])
  }

  if (i == 3) {
    text(x = -18, y = 1.45, labels = scenario[i])
  }

  if (i == 4) {
    text(x = -11, y = 1.45, labels = scenario[i])
  }

}

dev.off()
