library(copula)
library(doParallel)
library(doRNG) #for doParallel
library(ggplot2)


main <- function() {
  sink(file = "2019-02-21-simulation-log.txt")
  seed <- seed(525925) #525925
  tau <- c(0, 0.25, 0.5, 0.75) #Kendall's tau
  B <- c(10, 50, 100) #sample size of test statistics
  m <- c(6, 30, 60) #dimension
  default <- list(tau = 0.5, B = 100, m = 6) #default setting

  compare_generator(tau, B, m, default, seed)
}


#functions ####

calculate_phi.GNZ <- function(w, p, m) {
  n <- length(w)
  r_max <- 1000
  #r[1] is set to r_max => phi.GNZ is normalized with phi.GNZ(0) = r_max
  #this influences the behaviour of uniroot, but it does not change the estimated copula
  r <- c(r_max, numeric(n - 1)) #will be decreasing
  for (j in 2:n) {
    f <- function(x) sum((1 - x/r[1:(j - 1)])^(m - 1)*p[1:(j - 1)]) - w[j]
    #tol related errors are possible => try()
    r[j] <- uniroot(f, c(0, r[j - 1]), maxiter = 1000, tol = .Machine$double.eps)$root
  }
  r
}

calculate_phi.MC<- function(U, grid, n, m, M = 100) {
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
  phi_hat <- matrix(nrow = M, ncol = length(grid))
  for (b in 1:M) {
    U_ <- U[, sample(1:m, 2)] #submatrix
    w_ <- sort(apply(U_, 1, function(x) sum(colSums(t(U_) < x) == 2)))/(n + 1)
    w <- unique(w_) #depends on U_
    p <- vapply(w, function(x) mean(w_ == x), NA_real_)
    phi_hat_ <- phi(w, p)
    phi_hat[b, ] <- approx(c(w, 1), c(phi_hat_, 0), xout = grid)$y
  }
  phi_hat <- colMeans(phi_hat)
}

compare_generator <- function(tau, B, m, default, seed) {
  registerDoParallel(cores = max(detectCores() - 1, 1))
  results <- list()
  results <- perform_simulation(tau, default, seed, results)
  results <- perform_simulation(B, default, seed, results)
  results <- perform_simulation(m, default, seed, results)
  closeAllConnections()
  plot_results(results)
  results
}

#Performs a simulation regarding "var" and updates the list "results".
perform_simulation <- function(var, default, seed, results = list()) {
  set.seed(seed)
  varname <- deparse(substitute(var))
  results_update <- data.frame(matrix(, ncol = 5))
  for (var_ in var) {
    t <- Sys.time()
    print(paste(varname, " = ", var_, sep = ""))
    default[[varname]] <- var_

    results__ <- perform_test(default$tau, default$B, default$m)
    results_ <- data.frame(matrix(nrow = 3*length(results__$u), ncol = 5))
    GNZ <- 1:length(results__$u)
    MC <- (length(results__$u)+1):(2*length(results__$u))
    phi <- (2*length(results__$u)+1):(3*length(results__$u))

    results_[GNZ, 1] <- "GNZ"
    results_[MC, 1] <- "MC"
    results_[phi, 1] <- "phi"
    results_[, 2] <- var_

    results_[GNZ, 3] <- results__$mean$GNZ
    results_[MC, 3] <- results__$mean$MC
    results_[phi, 3] <- results__$mean$phi

    results_[GNZ, 4] <- results__$sd$GNZ
    results_[MC, 4] <- results__$sd$MC
    results_[phi, 4] <- results__$sd$phi

    results_[GNZ, 5] <- results__$u
    results_[MC, 5] <- results__$u
    results_[phi, 5] <- results__$u

    results_update <- rbind(results_update, results_)
    print(Sys.time() - t)
  }
  results_update <- results_update[-1, ]
  colnames(results_update) <- c("type", varname, "mean", "sd", "u")
  #update results
  results <- c(results, eval(parse(text = paste("list(results_", varname, " = results_update)",sep = ""))))
  save(results, file = "2019-02-21-simulation.RData")
  results
}

perform_test <- function(tau, B, m) {
  L <- 1000 #number of repetitions
  alpha <- 0.05
  n_try <- 2
  grid <- seq(0, 1, 1/(B + 1))
  results <- foreach(l = 1:L, .export = lsf.str(.GlobalEnv), .packages = "copula", .combine = rbind) %dorng% {
    result <- rep(NA, times = 2*length(grid))
    for (try_ in 1:n_try) {
      try(silent = TRUE, {
        C <- gumbelCopula(iTau(gumbelCopula(), tau), dim = m, use.indepC = "FALSE")
        U <- rCopula(B, C)
        w_ <- sort(apply(U, 1, function(x) sum(colSums(t(U) < x) == m)))/(B + 1)
        w <- unique(w_)
        p <- vapply(w, function(x) mean(w_ == x), NA_real_)
        phi_hat_ <- calculate_phi.GNZ(w, p, m)
        phi_hat.GNZ <- approx(c(w, 1), c(phi_hat_, 0), xout = grid)$y
        cop.GNZ <- approx(phi_hat.GNZ, grid, xout = m*phi_hat.GNZ)$y
        cop.GNZ[m*phi_hat.GNZ >= phi_hat.GNZ[1]] <- 0

        phi_hat.MC <- calculate_phi.MC(U, grid, B, m)
        cop.MC <- approx(phi_hat.MC, grid, xout = m*phi_hat.MC)$y
        cop.MC[m*phi_hat.MC >= phi_hat.MC[1]] <- 0

        result <- c(GNZ = cop.GNZ, MC = cop.MC)
        break
      })
    }
    result
  }
  if (sum(is.na(results)) > 0) {
    print("Warning: NAs produced.")
  }
  C <- gumbelCopula(iTau(gumbelCopula(), tau), dim = m, use.indepC = "FALSE")
  GNZ <- 1:length(grid)
  MC <- (length(grid)+1):(2*length(grid))
  phi <- (2*length(grid)+1):(3*length(grid))
  results_mean <- list(GNZ = colMeans(results[, GNZ], na.rm = TRUE),
                       MC = colMeans(results[, MC], na.rm = TRUE),
                       phi = pCopula(matrix(grid, nrow = length(grid), ncol = m), C))
  results_sd <- list(GNZ = apply(results[, GNZ], 2, sd, na.rm = TRUE),
                     MC = apply(results[, MC], 2, sd, na.rm = TRUE),
                     phi = rep(0, length(grid)))
  list(mean = results_mean, sd = results_sd, u = grid)
}

theme_set(theme_bw(base_size = 20))
theme_update(plot.title = element_text(hjust = 0.5))
update_geom_defaults("point", list(size = 2))
update_geom_defaults("line", list(size = 0.5))

plot_results <- function(results) {
  # methods <- c("GNZ", "MC")
  for (results_ in results) {
    if (nrow(results_) > 1) {
      varname <- names(results_)[2]
      plot_simulation(results_[, c("type", varname, "mean", "u")], varname,
                      title = "Mean")
      plot_simulation(results_[, c("type", varname, "sd", "u")], varname,
                      title = "SD")
    }
  }
}

plot_simulation <- function(data, varname, title) {
  short_title <- names(data)[3]
  if (varname == "B") {
    phi <- data[, "type"] == "phi" #the true copula does not change in B => only plot once
    gg <- ggplot() +
      labs(x = "u", y = "Copula(u)", color = varname, linetype = "") +
      ggtitle(title) +
      geom_line(data = data[phi, ], aes_string("u", short_title, group = "type")) +
      geom_line(data = data[!phi, ], aes_string("u", short_title, color = paste0("factor(", varname, ")"), group = paste0("interaction(", varname, ", type)"), linetype = "type"))
  } else {
    gg <- ggplot(data, aes_string("u", short_title, color = paste0("factor(", varname, ")"), group = paste0("interaction(", varname, ", type)"), linetype = "type")) +
      labs(x = "u", y = "Copula(u)", color = varname, linetype = "") +
      ggtitle(title) +
      geom_line()
  }
  print(gg)
  ggsave(paste0("2019-02-21-simulation-", varname, "-", short_title, ".png"), width = 6.4, height = 6.4)
}

#Prints a given random seed or creates a new one.
seed <- function(seed = NULL) {
  if (is.null(seed))
    seed <- as.numeric(paste(floor(runif(6, 0, 10)), collapse = ""))
  print(list(seed = seed))
  seed
}

####