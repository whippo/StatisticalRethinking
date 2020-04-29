#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Statistical Rethinking                                                         ##
# Data are current as of 2020-04-21                                              ##
# Data source: Statistical Rethinking Text                                       ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2020-04-21                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY: Problems and exercises from the book/lecutre series "Statistical
# Rethinking" by Richard McElreath (2nd edition). Sections are divided by
# chapter with annotation by Ross Whippo.


# Required Files (check that script is loading latest version):
# FILE.csv

# Associated Scripts:
# FILE.R

# TO DO 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                            ####
#                                                                                 +
# RECENT CHANGES TO SCRIPT                                                        +
# LOAD PACKAGES                                                                   +
# READ IN AND PREPARE DATA                                                        +
# MANIPULATE DATA                                                                 +
#                                                                                 +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RECENT CHANGES TO SCRIPT                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2020-04-21 Script created

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CHAPTER 2 - SMALL WORLDS AND LARGE WORLDS                                    ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 2.1 The garden of forking data ####


# 2.1.3 How to calculate the plausibility of an outcome based on prior information

# given 4 values: some combination of 0 and 1, what is the actual combination?
# you randomly draw (with replacment) three values and get: 1, 0, 1
# 0, 0, 0, 0 = 0 ways to get observation
# 1, 0, 0, 0 = 3 ways 
# 1, 1, 0, 0 = 8 ways
# 1, 1, 1, 0 = 9 ways
# 1, 1, 1, 1 = 0 ways

# You can update predictions of what the combination of values are by making a
# new draw of the same type of data. An additional draw produces: 1, now
# multiply the number of ways to draw 1's for each possible combination [W new]
# by prior observations [W prior].
# 0, 0, 0, 0 = 0 [D prior] * 0 [D new] = 0 
# 1, 0, 0, 0 = 3 [D prior] * 1 [D new] = 3
# 1, 1, 0, 0 = 8 [D prior] * 2 [D new] = 16
# 1, 1, 1, 0 = 9 [D prior] * 3 [D new] = 27
# 1, 1, 1, 1 = 0 [D prior] * 0 [D new] = 0

# You can also update predictions with diffrent types of information. Say that
# you are told that a set with a single '1' is 3x as likely as three '1's and
# that two '1's are twice as likely as three '1's. Again, you can multiply by
# the W prior.
# 0, 0, 0, 0 = 0  [D prior] * 0 [known likelihood] = 0 
# 1, 0, 0, 0 = 3  [D prior] * 3 [known likelihood] = 9
# 1, 1, 0, 0 = 16 [D prior] * 2 [known likelihood] = 32 
# 1, 1, 1, 0 = 27 [D prior] * 1 [known likelihood] = 27
# 1, 1, 1, 1 = 0  [D prior] * 0 [known likelihood] = 0

# 2.1.3 From counts to probablility

# To define probability of drawing a '1', you need to define the proportion of
# numbers that are '1' in each combination.
# 0, 0, 0, 0 = 0 
# 1, 0, 0, 0 = 0.25 
# 1, 1, 0, 0 = 0.50
# 1, 1, 1, 0 = 0.75
# 1, 1, 1, 1 = 1
# And given:
# D new = 1, 0, 1

# plausibility of p after D new ∝ ways p can produce D new × prior plausibility
# of p

# Standardize so that values are between 1 and 0 by dividing by sum of products.

# plausibility of p after D new = (ways p can produce D new × prior plausibility
# of p) / sum of products

# in R for all combinations:
ways <- c(0, 3, 8, 9, 0)
ways/sum(ways)
# [1] 0.00 0.15 0.40 0.45 0.00 
# these represent the posterior probabilities given the new information


# 2.3 Components of the model ####

# You want to know how much water is on the earth by tossing a globe into the
# air and counting the number of times your finger lands on water 'W' or land
# 'L'.

# 2.3.2.1 Observed variables

# If we say that all tosses are independent, and that the probability of landing
# on water 'W' is the same on every toss, it fits the binomial distribution.
# So for 9 tosses that resulted in 6 'W's and 3 'L's:
dbinom(6, size = 9, prob = 0.5)
# [1] 0.1640625 
# This is the relative number of ways to get 6 water given the conditions. Same
# way that posterior probabilities for 0 and 1 were calculated above.

# 2.4.3 Grid approximation

# In order to generate priors (when you can't calculate it easily like in the
# above examples), you must generate priors that are conditioned on the data.
# One way to do this for simple examples is by Grid approximation. It is a
# numerical technique for computing posterior distributions.

# Steps for grid approximation:
# 1) define the grid (how many points? list paramter values on grid)
# 2) compute value of the prior at each parameter value on the grid
# 3) compute the likelihood at each parameter value
# 4) compute the undstandardized posterior at each parameter value: prior *
# likelihood
# 5) standardize the posterior but dividing each value by the sum of all values
# (so that all values sum to 1)

# Example (using globe tossing example above):

# define grid
p_grid <- seq(from = 0, to = 1, length.out = 100)

# define prior
prior <- rep(1, 100) # a 'flat' prior, all equal probability
# prior <- ifelse(p_grid < 0.5, 0, 1) # a prior with no values less than 0.5
#prior <- exp(-5 * abs(p_grid - 0.5)) # a prior with a sharp peak at 0.5

# compute likelihood at each value in the grid
likelihood <- dbinom(6, size = 9, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# visualize the posterior
plot(p_grid, posterior, type = "b", xlab = "probability of water", ylab = "posterior probability")
mtext("100 points")

# 2.4.4 Quadratic approximation

# Grid approximation is bad for anything that involves multiple parameters
# because they quickly inflate and require more computing power. Quadratic
# approximation allows you to identify peaks in your distribution, then apply a
# Gaussian distribution to just the peak. These peaks are often nearly 'normal',
# so the distribution applies

# Example (using globe tossing example):

library(rethinking)
globe.qa <- quap(
  alist(
    W ~ dbinom(W + L, p), # binomial likelihood
    p ~ dunif(0, 1) # uniform prior
  ),
  data = list(W = 6, L = 3)
)

# display summary of quadratic approximation

precis(globe.qa)
#   mean   sd 5.5% 94.5%
# p 0.67 0.16 0.42  0.92

# analytical calculation to validate the above approximation
W <- 6
L <- 3
curve(dbeta(x, W + 1, L + 1), from = 0, to = 1)
#quadratic approximation
curve(dnorm(x, 0.67, 0.16), lty = 2, add = TRUE)

# with an N = 9, the approximation isn't great (assigns positive probability to
# 100% W which we know isn't true), but as sample size increases with same ratio
# of L:W, the curve improves until it's almost the same at N = 36

# 2.4.5 Markov chain Monte Carlo

# Useful for multilevel models that would take far to long to calculate a
# posterior any other way. MCMC draws samples from the posterior and the
# histograms of these frequencies correspond to the posterior plausibilities.

# Quick example how it works:

n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3
for (i in 2:n_samples) {
  p_new <- rnorm(1, p[i-1], 0.1)
  if(p_new < 0) p_new <- abs(p_new)
  if(p_new > 1) p_new <- 2 - p_new
  q0 <- dbinom(W, W+L, p[i-1])
  q1 <- dbinom(W, W+L, p_new)
  p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])
}

# values of p are samples from the posterior distribution 

# validate with algorithm

dens(p ,xlim=c(0,1))
curve(dbeta(x, W+1, L+1), lty=2, add=TRUE)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CHAPTER 3                                                                    ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# There are vampires in the world! And a blood test will tell you if you are one.

# Pr(positive test result|vampire) = 0.95; given that you are a vampire, the
# chance of testing positive is 95%

# Pr(positive test result|mortal) = 0.01; given that you are mortal, the chance
# of getting a false positive is 1%

# Pr(vampire) = 0.001; vampires are rare, only 0.1% of the population.

# Suppose someone tests positive for vampirism, what is the chance that they are
# actually a vampire? Must consider the true and false detection rates.

# Pr(vampire|positive) = Pr(positive|vampire) * Pr(vampire) / Pr(positive)
# Where Pr(positive) is the average probability of testing positive.
# Pr(positive) = Pr(positive|vampire) * Pr(vampire) + 
#                Pr(positive|mortal) * (1-Pr(vampire))

Pr_Positive_Vampire <- 0.95
Pr_Positive_Mortal <- 0.01
Pr_Vampire <- 0.001
Pr_Positive <- Pr_Positive_Vampire * Pr_Vampire +
               Pr_Positive_Mortal * (1-Pr_Vampire)
(Pr_Vampire_Positive <- Pr_Positive_Vampire * Pr_Vampire / Pr_Positive)
#[1] 0.08683729; only an 8.7% chance that they are really a vampire!

# Presented in other words, suppose you are told:
# 1) in a population of 100,000 people, 100 are vampires
# 2) of the 100 who are vampires, 95 of them will test positive 
# of the 99,900 mortals, 999 of them will test positive

# Now, what proportion tested will actually be vampires?

# total positive tests = 95 + 999 = 1094, so...

# Pr(vampire|positive) = 95/1094 = 0.08683729


# 3.1 Sampling from a grid-approximate posterior ####

# Use the globe tossing example above:
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, size = 9, prob = p_grid)
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

# If you draw sapmles from this posterior distribution, you will get samples in
# proportion to their probability

# Draw 10,000 samples:

samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

plot(samples)
library(rethinking)
dens(samples)


# 3.2 Sampling to summarize ####

# 3.2.1 Intervals of defined boundaries

# What is the posterior probability that the proportion of water is <0.5?

# add up posterior probability where p < 0.5
sum(posterior[p_grid < 0.5])
# [1] 0.1718746

# summing the posterior in this way only works for single parameters. For
# multiple parameters you can sum the sample, then divide by number of samples.

sum(samples < 0.5) / 1e4
#[1] 0.1737

# how much posterior probability lies between 0.5 and 0.75?

sum(samples > 0.5 & samples < 0.75) / 1e4
# [1] 0.6107


# 3.2.2 Invervals of defined mass

 # What is the lower bound of 80% of the data? (Compatibility interval [aka:
 # Confidence interval])
quantile(samples, 0.8)
# 0.7567568; 80% of the data lies below the value 0.76

# WHere is the middle 80% of data located?
quantile(samples, c(0.1, 0.9))
# 0.4504505 0.8119119; the middle 80% of data is between 0.45 and 0.81
# these are 'percentile intervals' PI, but they aren't great for skewed data

# what if you got water on three of three tosses?
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(3, size = 3, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)
plot(samples)
dens(samples)

# calculate the PI for central 50%
PI(samples, prob = 0.5)
# 0.7047047 0.9299299 
# it left out the most probable values near 1

# Instead, calculate where the 50% of highest density is with the Highest
# Posterior Density Interval (HPDI)

HPDI(samples, prob = 0.5)
# 0.8378378 0.9989990; the 50% highest density of probability lies between 84%
# and 99%








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# HOMEWORK - WEEK 1                                                            ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. Construct a posterior distribution for the globe exercies as if there had
# been 8 water in 15 tosses using grid approximation. Use the same flat prior.

# define grid
p_grid <- seq(from = 0, to = 1, length.out = 100)

# define prior
prior <- rep(1, 100) # a 'flat' prior, all equal probability

# compute likelihood at each value in the grid
likelihood <- dbinom(8, size = 15, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# visualize the posterior
plot(p_grid, posterior, type = "b", xlab = "probability of water", ylab = "posterior probability")
mtext("100 points")

# 2. Do again, but this time use a prior that is zero below p = 0.5 and a
# constant above p = 0.5, corresponding to the knowledge that the earth is
# mostly water. What difference does a better prior make?

# define grid
p_grid <- seq(from = 0, to = 1, length.out = 100)

# define prior
prior <- ifelse(p_grid < 0.5, 0, 1) # a prior with no values less than 0.5

# compute likelihood at each value in the grid
likelihood <- dbinom(8, size = 15, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# visualize the posterior
plot(p_grid, posterior, type = "b", xlab = "probability of water", ylab = "posterior probability")
mtext("100 points")

# 3. How would you estimate the earth's proportion of water precisely, obtaining
# the 99th percentile interval of the posterior distribution of p to be only
# 0.05 wide? This means the distance between upper and lower bound should only
# be 0.05. How many times will you have to toss the globe?

############### SUBSECTION HERE

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####