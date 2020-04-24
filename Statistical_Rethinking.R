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
# 100% W    which we know isn't true), but as sample size increases with same ratio
# of L:W, the curve improves until it's almost the same at N = 36


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


############### SUBSECTION HERE

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####