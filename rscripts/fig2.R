#AJH modified version of RZF code for example in Rannala, 2002

library(RColorBrewer)
palette(brewer.pal(4, "Set2"))
cols <- brewer.pal(4, "Set2")

#number of simulated data points
n =  50

#true value of parameter1
a_true = 0.5

#true value of parameter2
b_true = 0.25

#sum of true values
diff_true <- a_true - b_true

#simulate data with gamma distribution
#Density, distribution function, quantile function and random generation
#for the Gamma distribution with parameters shape and scale.

simulated_data <- rgamma(n, #number of simulated data points
                         shape = 2, #shape of gamma distribution
                         rate = diff_true) #scale


#look at simulated data
plot(simulated_data)

#rate for exponential distribution
alpha = 0.5

#draw 2 values from exponential disribution
m1 <- rexp(1, 1 / alpha)
m2 <- rexp(1, 1 / alpha)

#function to calculate likelihood
likelihood <- function(lambda, mu, datos) {
  if ((lambda - mu) < 0) {
    like <- 0
    return(like) #fail
  } else{
    like <-
      prod(dgamma(datos, 2, rate = (lambda - mu))) #product of all vectors of elements in gamma dist
    return(like)
  }
}

#function to calculate prior
prior <- function(lambda, mu) {
  if ((lambda - mu) < 0) {
    prior.val <- 0
    return(prior.val)
  } else{
    prior.val <- dexp(lambda, 1 / alpha) * dexp(mu, 1 / alpha)
    return(prior.val)
  }
}

#number of generations to run chain
generations <- 5000

#limits on sampling (+ or - this value)
delta <- 0.25

#prepare output matrix
output <- matrix(rep(0, 6 * generations), ncol = 6)


for (i in 1:generations) {
  #modify params (step)
  lambda_prime <- m1 + runif(1, -delta, delta)
  mu_prime <- m2 + runif(1, -delta, delta)

  #calculate ratio of likelihoods of new_vals/old_vals
  like_odds <-
    likelihood(lambda_prime, mu_prime, simulated_data) / likelihood(m1, m2, simulated_data)

  #calculate ratio of prior of new_vals/old_vals
  prior_odds <- prior(lambda_prime, mu_prime) / prior(m1, m2)

  #calculate posterior odds
  R <- like_odds * prior_odds

  #randomly draw from uniform dist
  u <- runif(1)

  #if posterior odds are greater than random value, keep new values of parameter
  if (u < R) {
    m1 = lambda_prime
    m2 = mu_prime
  }

  #calculate posterior (likelihood*prior)
  posterior <- likelihood(m1, m2, simulated_data) * prior(m1, m2)

  #store output
  output[i, 1] <- i
  output[i, 2] <- prior(m1, m2)
  output[i, 3] <- likelihood(m1, m2, simulated_data)
  output[i, 4] <- posterior
  output[i, 5] <- m1
  output[i, 6] <- m2
}

#format output data
output <- data.frame(output)
names(output) <-
  c("iteration", "prior", "likelihood", "posterior", "lambda", "mu")

pdf(
  "figures/fig2a-c.pdf",
  width = 5,
  height = 9
)
par(mar = c(4, 4, 2, 2))
#plot values of parameters in chain
par(mfrow = c(3, 1))
plot(
  output$lambda ~ output$iteration,
  type = 'l',
  col = 4,
  ylim = c(0, 1.5),
  xlab = "Generation",
  ylab = "lambda",
  bty = 'l'
)
abline(h = a_true, lty = 2)
mtext(
  "a)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

plot(
  output$mu ~ output$iteration,
  type = 'l',
  col = 2,
  ylim = c(0, 1.5),
  xlab = "Generation",
  ylab = "mu",
  bty = 'l'
)
abline(h = b_true, lty = 2)
mtext(
  "b)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

plot((output$lambda - output$mu) ~ output$iteration,
     type = 'l',
     col = 3,
     ylim = c(0, 1.5),
     xlab = "Generation",
     ylab = "lambda - mu",
     bty = 'l'
)
abline(h = a_true - b_true, lty = 2)
mtext(
  "c)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

dev.off()

sd(output$prior)

#uncomment to save to pdf
#pdf(
#  "figures/fig2_prior_lik_posterior.pdf",
#  width = 5,
#  height = 9
#)

#plot prior, likelihood and posterior
par(mfrow = c(3, 1))

par(mar = c(4, 4, 2, 2))

plot(
  output$iteration,
  output$prior,
  col = 1,
  type = 'l',
  xlab = "Generation",
  ylab = "Prior",
  ylim = c(
    min(output$prior) - sd(output$prior),
    max(output$prior) + sd(output$prior)
  )
)
mtext(
  "a)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

plot(
  output$iteration,
  output$likelihood,
  col = 2,
  type = 'l',
  xlab = "Generation",
  ylab = "Likelihood",
  ylim = c(
    min(output$likelihood) - sd(output$likelihood),
    max(output$likelihood) + sd(output$likelihood)
  )
)
mtext(
  "b)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

plot(
  output$iteration,
  output$posterior,
  col = 3,
  type = 'l',
  xlab = "Generation",
  ylab = "Posterior",
  ylim = c(
    min(output$posterior) - sd(output$posterior),
    max(output$posterior) + sd(output$posterior)
  )
)
mtext(
  "c)",
  side = 3,
  line = 1,
  adj = -0.1,
  cex = 0.8
)

#uncomment to save to pdf
#dev.off()

par(mfrow = c(1, 1))

#vector across range of values interested in
values_m <- seq(0, 1.2, 0.005)

#Create a data frame from all combinations of the supplied vector
contour_mat <- expand.grid(values_m, values_m)

#prepare for storing likelhood
likelihood_vals <- rep(0, dim(contour_mat)[1])

#for each pair of values calculate likelihood
for (i in 1:dim(contour_mat)[1]) {
  likelihood_vals[i] <-
    likelihood(contour_mat[i, 1], contour_mat[i, 2], simulated_data)
}

#get location of maximum likelihood value
maximum <- which.max(likelihood_vals)

#convert raw likelihoods to relative of ML
relative = likelihood_vals / likelihood_vals[maximum]
contour_mat <- cbind(contour_mat, likelihood_vals)
contour_mat <- cbind(contour_mat, relative)

library("ggplot2")
jpeg(
  "figures/fig2e.jpeg",
  width = 6,
  height = 6,
  units = "in",
  res = 500
)

ggplot(contour_mat, aes(x = Var1, y = Var2, z = relative)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = relative)) +
  stat_contour(bins = 10,
               color = cols[1],
               size = 0.5)  +
  scale_fill_gradient(
    low = "white",
    high = cols[1],
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  xlab("lambda") +
  ylab("mu") +
  xlim(0.2, 1) +
  ylim(0, 0.8) +
  guides(fill = guide_colorbar(title = "Relative Likelihood")) +
  geom_point(
    aes(x = a_true, y = b_true),
    colour = cols[2],
    pch = 16,
    size = 5,
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey"),
    legend.position = c(0.2, 0.8)
  )

dev.off()

pdf(
  "figures/fig2d.pdf",
  width = 6,
  height = 6
)
ggplot(output, aes(x = lambda, y = mu)) + geom_point(shape = 19, color = cols[1]) +
  geom_smooth(
    method = lm,
    se = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  xlim(0.2, 1) +
  ylim(0, 0.8) +
  geom_point(
    aes(x = a_true, y = b_true),
    colour = cols[2],
    pch = 16,
    size = 5
  ) +
  theme(
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(colour = "grey"),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )
dev.off()
