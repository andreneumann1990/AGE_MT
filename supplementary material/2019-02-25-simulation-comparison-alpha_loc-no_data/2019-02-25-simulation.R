library(copula)
library(doParallel)
library(doRNG) #for doParallel
library(ggplot2)

main <- function() {
  sink(file = "2019-02-25-simulation-log.txt")
  seed <- seed(717801) #717801
  tau <- 0:19/20 #Kendall's tau
  B <- c(10*1:10, 100*2:10) #size of the test statistics sample under LFC
  m <- c(2:9, 10*1:9, 100*1:10) #number of hypotheses
  pi0 <- 1:5/6 #proportion of true null hypotheses
  effect <- 1:4/10 #effect size under alternatives
  default <- list(tau = 0.5, B = 100, m = 6, pi0 = 0.5, effect = 0.2) #default setting

  compare_generator(tau, B, m, pi0, effect, default, seed)
}


#functions ####

calculate_phi.GNZ <- function(w, p, d) {
  n <- length(w)
  r_max <- 1000
  #r[1] is set to r_max => phi.GNZ is normalized with phi.GNZ(0) = r_max
  #this influences the behaviour of uniroot, but it does not change the estimated copula
  r <- c(r_max, numeric(n - 1)) #will be decreasing
  for (j in 2:n) {
    f <- function(x) sum((1 - x/r[1:(j - 1)])^(d - 1)*p[1:(j - 1)]) - w[j]
    #tol related errors are possible => try()
    r[j] <- uniroot(f, c(0, r[j - 1]), maxiter = 1000, tol = .Machine$double.eps)$root
  }
  r
}

calculate_phi.MC<- function(U, possible_w, B, m, MC = 100) {
  phi <- function(w, p) {
    n <- length(w)
    K_hat <- numeric(n)
    K_hat[1] <- p[1]
    for (i in 2:n) {
      K_hat[i] <- K_hat[i - 1] + p[i]
    }
    r_max <- 1000
    r <- c(r_max, numeric(n - 1)) #will be decreasing
    for (j in 2:n) {
      #this is essentially the same what calculate_phi.GNZ does in the case m = 2 but directly without using uniroot
      r[j] <- r[j - 1] * (w[j] - K_hat[j - 1])/(w[j - 1] - K_hat[j - 1])
    }
    r
  }
  phi_hat <- matrix(nrow = MC, ncol = B + 1)
  for (b in 1:MC) {
    U_ <- U[, sample(1:m, 2)] #submatrix
    w_ <- sort(apply(U_, 1, function(x) sum(colSums(t(U_) < x) == 2)))/(B + 1)
    w <- unique(w_) #depends on U_
    p <- vapply(w, function(x) mean(w_ == x), NA_real_)
    phi_hat_ <- phi(w, p)
    phi_hat[b, ] <- approx(c(w, 1), c(phi_hat_, 0), xout = possible_w)$y
  }
  phi_hat <- colMeans(phi_hat)
}

compare_generator <- function(tau, B, m, pi0, effect, default, seed) {
  registerDoParallel(cores = max(detectCores() - 1, 1))
  results <- list()
  results <- perform_simulation(tau, default, seed, results)
  results <- perform_simulation(B, default, seed, results)
  results <- perform_simulation(m, default, seed, results)
  results <- perform_simulation(pi0, default, seed, results)
  results <- perform_simulation(effect, default, seed, results)
  closeAllConnections()
  plot_results(results)
  results
}

#Performs a simulation regarding "var" and updates the list "results".
perform_simulation <- function(var, default, seed, results = list()) {
  set.seed(seed)
  varname <- deparse(substitute(var))
  l <- length(var)
  types <- c("GNZ", "MC", "phi", "Bonf")
  results_ <- data.frame(matrix(nrow = length(types)*l, ncol = 4))
  names(results_) <- c("type", varname, "mean", "sd")
  k <- 0
  for (type in types) {
    results_[(k*l+1):((k+1)*l), 1] <- type
    k <- k+1
  }
  results_[, 2] <- rep(var, times = length(types))
  i <- 1
  for (var_ in var) {
    t <- Sys.time()
    print(paste(varname, " = ", var_, sep = ""))
    default[[varname]] <- var_
    results__ <- perform_test(default$tau, default$B, default$m, default$pi0, default$effect)
    for (k in 1:length(types)) {
      results_[i+(k-1)*l, 3:4] <- results__[c(k, length(types)+k)]
    }
    i <- i + 1
    print(Sys.time() - t)
  }
  #update results
  results <- c(results, eval(parse(text = paste("list(results_", varname, " = results_)",sep = ""))))
  save(results, file = "2019-02-25-simulation.RData")
  results
}

perform_test <- function(tau, B, m, pi0, effect) {
  L <- 1000 #number of repetitions
  alpha <- 0.05
  n_try <- 2
  m0 <- ceiling(m*pi0)
  possible_w <- c(0:(B - 1)/(B + 1), 1)
  results <- foreach(l = 1:L, .export = lsf.str(.GlobalEnv), .packages = "copula", .combine = rbind) %dorng% {
    result <- rep(NA, 4)
    for (try_ in 1:n_try) {
      try(silent = TRUE, {
        C <- gumbelCopula(iTau(gumbelCopula(), tau) + (m-m0)/m*effect, dim = m, use.indepC = "FALSE")
        T <- qnorm(rCopula(1, C))
        T[, (m0+1):m] <- T[, (m0+1):m] + effect
        C_star <- gumbelCopula(iTau(gumbelCopula(), tau), dim = m, use.indepC = "FALSE")
        T_star <- qnorm(rCopula(B, C_star)) # sample of T under LFC
        w_ <- sort(apply(T_star, 1, function(x) sum(colSums(t(T_star) < x) == m)))/(B + 1)
        w <- unique(w_)
        p <- vapply(w, function(x) mean(w_ == x), NA_real_)
        phi_hat_ <- calculate_phi.GNZ(w, p, m)
        phi_hat.GNZ <- approx(c(w, 1), c(phi_hat_, 0), xout = possible_w)$y[-1]
        phi_hat.GNZ <- phi_hat.GNZ/phi_hat.GNZ[1]
        phi_of_alpha_loc.GNZ <- 1/m*approx(possible_w[-1], phi_hat.GNZ, xout = 1 - alpha)$y
        alpha_loc.GNZ <- 1 - approx(phi_hat.GNZ, possible_w[-1], xout = phi_of_alpha_loc.GNZ)$y

        phi_hat.MC <- calculate_phi.MC(T_star, possible_w, B, m)[-1]
        phi_hat.MC <- phi_hat.MC/phi_hat.MC[1]
        phi_of_alpha_loc.MC <- 1/m*approx(possible_w[-1], phi_hat.MC, xout = 1 - alpha)$y
        alpha_loc.MC <- 1 - approx(phi_hat.MC, possible_w[-1], xout = phi_of_alpha_loc.MC)$y

        alpha_loc.phi <- 1 - psi(C, 1/m*iPsi(C, 1 - alpha))

        result <- c(GNZ = alpha_loc.GNZ, MC = alpha_loc.MC, phi = alpha_loc.phi, Bonf = alpha/m)
        break
      })
    }
    result
  }
  if (sum(is.na(results)) > 0) {
    print("Warning: NAs produced.")
  }
  c(mean = apply(results, 2, mean), sd = apply(results, 2, sd))
}

theme_set(theme_bw(base_size = 20))
theme_update(plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size = 2))
update_geom_defaults("line", list(size = 1.2))

plot_results <- function(results) {
  for (results_ in results) {
    if (nrow(results_) > 1) {
      varname <- names(results_)[2]
      plot_simulation(results_[, c("type", varname, "mean")], varname,
                      title = "Mean alpha_loc")
      plot_simulation(results_[, c("type", varname, "sd")], varname,
                      title = "SD alpha_loc")
    }
  }
}

plot_simulation <- function(data, varname, title) {
  short_title <- names(data)[3]
  gg <- ggplot(data, aes(data[, varname], data[, short_title], color = type)) +
    labs(x = parse(text = varname), y = "", colour = "") +
    ggtitle(title) +
    geom_line()
  ggsave(paste("2019-02-25-simulation-", varname, "-", short_title,
               ".png", sep = ""), width = 6.4, height = 6.4)
  print(gg)
}

#Prints a given random seed or creates a new one.
seed <- function(seed = NULL) {
  if (is.null(seed))
    seed <- as.numeric(paste(floor(runif(6, 0, 10)), collapse = ""))
  print(list(seed = seed))
  seed
}

####