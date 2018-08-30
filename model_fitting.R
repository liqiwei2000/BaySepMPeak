# ***README***
# The following script is used to fit MeRIP-Seq data (in multiple conditions) to the 
# Bayesian zero-inflated negative binomial linear regression model with spatial variable
# selection proposed in the submitted manuscript titled "..."

# Before running the following code, please first load MeRIP-Seq data. The necessary 
# inputs should be 
# (1) a n-by-W count matrix Y, where n is the number of samples and W is the number 
# of bins
# (2) a design matrix X, indiating the allocation of IP/case samples.
# (3) a n-dimensional vector, indicating the size factor for each sample.

# Note that the notations in the code and data follow the notations in the manuscript.
# ***END***

# Load functions
Rcpp::sourceCpp('core_2s.cpp');



# ========================================================================================
# ========================================================================================
# Load data
# ========================================================================================
# ========================================================================================
Y <- read.table("real_data_example.txt", header = TRUE);
Y <- as.matrix(t(Y));
X <- read.table("design_matrix.txt", header = TRUE);
X <- as.matrix(X);
s <- read.table("size_factors.txt", header = FALSE, stringsAsFactors = FALSE)[, 2];
n <- dim(Y)[1];
W <- dim(Y)[2];
K <- (dim(X)[2] - 1)/2;



# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
iter <- 100*W; 



# ========================================================================================
# ========================================================================================
# Load hyperparameters, please adjust the other hyperparameters in "zinb_hmm_x2.cpp"
# ========================================================================================
# ========================================================================================
# Hyperparameters for Markov random field prior
e <- -4.6;
f <- 2.3;
# Hyperparameters for Inverse-gamma prior
a <- 2;
b <- 1;

# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm
# ========================================================================================
# ========================================================================================
start_time <- proc.time();
M <- BaySeqMPeak(Y, X, s, e, f, a, b, iter);
end_time <- proc.time();
time <- end_time - start_time;
# The MCMC outputs are stored in M
# $flag:          the indicator of valid windows, TRUE means the total counts of the 
#                 corresponding window is less than 10 and the total number of non-zero 
#                 counts in the window is less than 2
# $pi:            the mean of MCMC samples of pi
# $H:             the marginal posterior probability of inclusion (PPI) of H
# $phi:           the mean of MCMC samples of phi
# $alpha_0_input: the mean of MCMC samples of alpha_0_input
# $alpha_0_ip:    the mean of MCMC samples of alpha_0_ip
# $alpha:         the mean of MCMC samples of alpha
# $gamma:         the marginal posterior probability of inclusion (PPI) of gamma
# $Theta:         the mean of MCMC samples of Theta
# $Xi:            the marginal posterior probability of inclusion (PPI) of Xi
# $Beta:          the mean of MCMC samples of Theta
# $Delta:         the marginal posterior probability of inclusion (PPI) of Delta
# $accept:        the acceptance rate for updating phi, alpha_0_input, alpha_0_ip,
#                 (alpha, Beta, Theta), and (gamma, Delta, Xi)



# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
par(mfrow = c(2, 2));
cutoff <- 0.5;
max <- max(Y/s, exp(M$alpha_0_ip + M$alpha), exp(M$alpha_0_ip + M$alpha + c(M$Beta)), exp(M$alpha_0_ip + M$alpha + c(M$Beta) + c(M$Theta)));

plot(NULL, ylim = c(0, max), xlim = c(1, W), xlab = "Bin", ylab = "Abundance", main = "Control input samples (3)");
lines(exp(M$alpha_0_input), col = 1);
index <- which(X[, 1] == 0 & X[, 2] == 0)
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1);
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff));
}

# max <- max(Y[index,]/s[index], exp(M$alpha_0_input + c(M$Beta)));
plot(NULL, ylim = c(0, max), xlim = c(1, W), xlab = "Bin", ylab = "Abundance", main = "Experimental input samples (3)");
lines(exp(M$alpha_0_input), col = 1);
lines(exp(M$alpha_0_input + c(M$Beta)), col = 4);
index <- which(X[, 1] == 0 & X[, 2] == 0);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
}
index <- which(X[, 1] == 0 & X[, 2] == 1);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = 4);
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = 4);
}
temp <- c(M$Delta);
index_temp <- which(temp >= cutoff);
if (length(index_temp) > 1) {
  break_temp <- rep(0, length(index_temp));
  count <- 1;
  break_temp[1] = count;
  for (i in 2:length(index_temp)) {
    if (index_temp[i - 1] + 1 != index_temp[i]) {
      count <- count + 1;
      break_temp[i] <- count;
    } else {
      break_temp[i] <- count;
    }
  }
  for (i in 1:max(break_temp)) {
    lines(index_temp[which(break_temp == i)], rep(max, sum(break_temp == i)), lwd = 5, col = 4);
  }
}

# max <- max(Y[index,]/s[index], exp(M$alpha_0_ip + M$alpha));
plot(NULL, ylim = c(0, max), xlim = c(1, W), xlab = "Bin", ylab = "Abundance", main = "Control IP samples (3)");
lines(exp(M$alpha_0_input), col = 1);
lines(exp(M$alpha_0_ip + M$alpha), col = 2);
index <- which(X[, 1] == 0 & X[, 2] == 0);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
}
index <- which(X[, 1] == 1 & X[, 2] == 0);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = 2);
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = 2);
}
temp <- M$gamma;
index_temp <- which(temp >= cutoff);
if (length(index_temp) > 1) {
  break_temp <- rep(0, length(index_temp));
  count <- 1;
  break_temp[1] = count;
  for (i in 2:length(index_temp)) {
    if (index_temp[i - 1] + 1 != index_temp[i]) {
      count <- count + 1;
      break_temp[i] <- count;
    } else {
      break_temp[i] <- count;
    }
  }
  for (i in 1:max(break_temp)) {
    lines(index_temp[which(break_temp == i)], rep(max, sum(break_temp == i)), lwd = 5, col = 2);
  }
}

# max <- max(Y[index,]/s[index], exp(M$alpha_0_ip + M$alpha + c(M$Beta) + c(M$Theta)));
plot(NULL, ylim = c(0, max), xlim = c(1, W), xlab = "Bin", ylab = "Abundance", main = "Experimental IP samples (3)");
lines(exp(M$alpha_0_input + c(M$Beta)), col = 4);
lines(exp(M$alpha_0_ip + M$alpha + c(M$Beta)), col = 2);
lines(exp(M$alpha_0_ip + M$alpha + c(M$Beta) + c(M$Theta)), col = 6);
index <- which(X[, 1] == 0 & X[, 2] == 0);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1));
}
index <- which(X[, 1] == 0 & X[, 2] == 1);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1));
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1));
}
index <- which(X[, 1] == 1 & X[, 2] == 0);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1));
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1));
}
index <- which(X[, 1] == 1 & X[, 2] == 1);
for (i in index) {
  points(1:W, Y[i,]/s[i], pch = 1, col = 6);
  points(1:W, Y[i,]/s[i], pch = 1 + 3*(M$H[i,] >= cutoff), col = 6);
}
temp <- c(M$Xi);
index_temp <- which(temp >= cutoff);
if (length(index_temp) > 1) {
  break_temp <- rep(0, length(index_temp));
  count <- 1;
  break_temp[1] = count;
  for (i in 2:length(index_temp)) {
    if (index_temp[i - 1] + 1 != index_temp[i]) {
      count <- count + 1;
      break_temp[i] <- count;
    } else {
      break_temp[i] <- count;
    }
  }
  for (i in 1:max(break_temp)) {
    lines(index_temp[which(break_temp == i)], rep(max, sum(break_temp == i)), lwd = 5, col = 6);
  }
}
par(mfrow = c(1, 1));