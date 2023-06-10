#' Pointwise sample mean functions
#'
#' The function \code{pointwise_sample_mean_fun()} calculates and draws the pointwise sample mean functions.
#'
#' @param x a list of length \eqn{\ell} with elements being \eqn{n\times p} matrices of data
#' corresponding to \eqn{n} functional observations measured in \eqn{p} design time points under given
#' experimental conditions.
#' @param plot a logical indicating of whether to draw the values of the pointwise sample mean functions.
#' The default is \code{TRUE}.
#' @param values a logical indicating of whether to return the values of the pointwise sample mean functions.
#' The default is \code{FALSE}.
#' @param type 1-character string giving the type of plot desired, the same as in the \code{matplot()} function.
#' The default is \code{"l"} for lines.
#' @param lty vector of line types, the same as in the \code{matplot()} function. The default is 1 (solid lines).
#' @param main a main title for the plot, the same as in the \code{plot()} function. The default is \code{Sample mean functions}.
#' @param ... other graphical parameters, the same as in the \code{matplot()} function.
#'
#' @return If \code{values = TRUE}, a matrix of values of the pointwise sample mean functions.
#'
#' @references Kurylo K., Smaga L. (2023) Functional repeated measures analysis of variance and its application. Preprint https://arxiv.org/abs/2306.03883
#'
#' @examples
#' # preparation of the DTI data set, for details see Kurylo and Smaga (2023)
#' library(refund)
#' data(DTI)
#' # MS patients
#' DTI_ms <- DTI[DTI$case == 1, ]
#' miss_data <- c()
#' for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
#' DTI_ms <- DTI_ms[-miss_data, ]
#' DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
#' xx <- vector("list", 4)
#' for (i in 1:4) {
#'   xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
#' }
#' xx[[1]] <- xx[[1]][-14, ]
#' xx[[3]] <- xx[[3]][-14, ]
#' yy <- xx
#' for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
#' # sample mean functions
#' oldpar <- par(mfrow = c(1, 1), mar = c(4, 4, 2, 0.1))
#' pointwise_sample_mean_fun(yy, values = FALSE,
#'                           col = 1:4, xlab = "t", ylab = "FA", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' legend(x = 36, y = 0.64, legend = 1:4, lty = 1, col = 1:4, title = "Visit")
#' par(oldpar)
#'
#' @export
pointwise_sample_mean_fun <- function(x, plot = TRUE, values = FALSE, type = "l", lty = 1,
                                      main = "Sample mean functions", ...) {
  p <- ncol(x[[1]])
  l <- length(x)
  means <- matrix(NA, nrow = p, ncol = l)
  for (i in seq_len(l)) {
    means[, i] <- colMeans(x[[i]])
  }
  if (plot) {
    matplot(means, type = type, lty = lty, main = main, ...)
  }
  if (values) {
    return(means)
  }
}

cn_value <- function(x) {
  n <- nrow(x[[1]])
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  return(n * sum((means_gr - means_all)^2))
}

dn_en_value <- function(x) {
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  k <- length(x)
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  means_sub <- matrix(0, nrow = n, ncol = p)
  for (ii in seq_len(k)) {
    means_sub <- means_sub + x[[ii]]
  }
  means_sub <- t(means_sub) / k
  SSA <- rowSums(n * (means_gr - means_all)^2)
  SSS <- rowSums(k * (means_sub - means_all)^2)
  SST <- rowSums(sapply(lapply(x, function(x) (t(x) - means_all)^2), rowSums))
  SSE <- SST - SSA - SSS
  f_point <- (SSA / (k - 1)) / (SSE / ((n - 1) * (k - 1)))
  f_point <- f_point[is.finite(f_point)]
  f_point <- ifelse(f_point < .Machine$double.eps, 0, f_point)
  return(c(sum(f_point), max(f_point)))
}

cn_value_boot <- function(x, means_gr_raw, means_all_raw) {
  n <- nrow(x[[1]])
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  return(n * sum((means_gr - means_gr_raw - means_all + means_all_raw)^2))
}

dn_en_value_boot <- function(x, means_gr_raw, means_all_raw) {
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  k <- length(x)
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  means_sub <- matrix(0, nrow = n, ncol = p)
  for (ii in seq_len(k)) {
    means_sub <- means_sub + x[[ii]]
  }
  means_sub <- t(means_sub) / k
  SSA_boot <- rowSums(n * (means_gr - means_gr_raw - means_all + means_all_raw)^2)
  SSA <- rowSums(n * (means_gr - means_all)^2)
  SSS <- rowSums(k * (means_sub - means_all)^2)
  SST <- rowSums(sapply(lapply(x, function(x) (t(x) - means_all)^2), rowSums))
  SSE <- SST - SSA - SSS
  f_point <- (SSA_boot / (k - 1)) / (SSE / ((n - 1) * (k - 1)))
  f_point <- f_point[is.finite(f_point)]
  f_point <- ifelse(f_point < .Machine$double.eps, 0, f_point)
  return(c(sum(f_point), max(f_point)))
}

#' Pointwise SSA test statistic
#'
#' The function \code{pointwise_ssa_test_statistic()} calculates and draws the pointwise SSA test statistic.
#'
#' @param x a list of length \eqn{\ell} with elements being \eqn{n\times p} matrices of data
#' corresponding to \eqn{n} functional observations measured in \eqn{p} design time points under given
#' experimental conditions.
#' @param plot a logical indicating of whether to draw the values of the pointwise SSA test statistic.
#' The default is \code{TRUE}.
#' @param values a logical indicating of whether to return the values of the pointwise SSA test statistic.
#' The default is \code{FALSE}.
#' @param type 1-character string giving the type of plot desired, the same as in the \code{plot()} function.
#' The default is \code{"l"} for lines.
#' @param ylab a label for the \eqn{y}-axis, the same as in the \code{plot()} function.
#' The default is the empty sign.
#' @param main a main title for the plot, the same as in the \code{plot()} function. The default is \code{SSA(t)}.
#' @param ... other graphical parameters, the same as in the \code{plot()} function.
#'
#' @details For details, see the documentation of the \code{rmfanova()} function or
#' the paper Kurylo and Smaga (2023).
#'
#' @return If \code{values = TRUE}, a vector of values of the pointwise SSA test statistic.
#'
#' @references Martinez-Camblor P., Corral N. (2011) Repeated Measures Analysis for Functional Data.
#' Computational Statistics & Data Analysis 55, 3244–3256.
#'
#' Kurylo K., Smaga L. (2023) Functional repeated measures analysis of variance and its application.
#' Preprint https://arxiv.org/abs/2306.03883
#'
#' @examples
#' # preparation of the DTI data set, for details see Kurylo and Smaga (2023)
#' library(refund)
#' data(DTI)
#' # MS patients
#' DTI_ms <- DTI[DTI$case == 1, ]
#' miss_data <- c()
#' for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
#' DTI_ms <- DTI_ms[-miss_data, ]
#' DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
#' xx <- vector("list", 4)
#' for (i in 1:4) {
#'   xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
#' }
#' xx[[1]] <- xx[[1]][-14, ]
#' xx[[3]] <- xx[[3]][-14, ]
#' yy <- xx
#' for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
#' # pointwise SSA test statistic
#' pointwise_ssa_test_statistic(yy, xlab = "t", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#'
#' @export
pointwise_ssa_test_statistic <- function(x, plot = TRUE, values = FALSE,
                                         type = "l", ylab = "", main = "SSA(t)", ...) {
  n <- nrow(x[[1]])
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  ssa <- n * rowSums((means_gr - means_all)^2)
  if (plot) {
    plot(ssa, type = type, ylab = ylab, main = main, ...)
  }
  if (values) {
    return(ssa)
  }
}

#' Pointwise F-type test statistic
#'
#' The function \code{pointwise_f_test_statistic()} calculates and draws the pointwise F-type test statistic.
#'
#' @param x a list of length \eqn{\ell} with elements being \eqn{n\times p} matrices of data
#' corresponding to \eqn{n} functional observations measured in \eqn{p} design time points under given
#' experimental conditions.
#' @param plot a logical indicating of whether to draw the values of the pointwise F-type test statistic.
#' The default is \code{TRUE}.
#' @param values a logical indicating of whether to return the values of the pointwise F-type test statistic.
#' The default is \code{FALSE}.
#' @param type 1-character string giving the type of plot desired, the same as in the \code{plot()} function.
#' The default is \code{"l"} for lines.
#' @param ylab a label for the \eqn{y}-axis, the same as in the \code{plot()} function.
#' The default is the empty sign.
#' @param main a main title for the plot, the same as in the \code{plot()} function. The default is \code{F(t)}.
#' @param ... other graphical parameters, the same as in the \code{plot()} function.
#'
#' @details For details, see the documentation of the \code{rmfanova()} function or
#' the paper Kurylo and Smaga (2023).
#'
#' @return If \code{values = TRUE}, a vector of values of the pointwise F-type test statistic.
#'
#' @references Kurylo K., Smaga L. (2023) Functional repeated measures analysis of variance and its application. Preprint https://arxiv.org/abs/2306.03883
#'
#' @examples
#' # preparation of the DTI data set, for details see Kurylo and Smaga (2023)
#' library(refund)
#' data(DTI)
#' # MS patients
#' DTI_ms <- DTI[DTI$case == 1, ]
#' miss_data <- c()
#' for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
#' DTI_ms <- DTI_ms[-miss_data, ]
#' DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
#' xx <- vector("list", 4)
#' for (i in 1:4) {
#'   xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
#' }
#' xx[[1]] <- xx[[1]][-14, ]
#' xx[[3]] <- xx[[3]][-14, ]
#' yy <- xx
#' for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
#' # pointwise F-type test statistic
#' pointwise_f_test_statistic(yy, xlab = "t", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#'
#' @export
pointwise_f_test_statistic <- function(x, plot = TRUE, values = FALSE,
                                       type = "l", ylab = "", main = "F(t)", ...) {
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  k <- length(x)
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  means_sub <- matrix(0, nrow = p, ncol = n)
  for (ii_s in seq_len(n)) {
    means_sub_temp <- 0
    for (jj_s in seq_len(k)) {
      means_sub_temp <- means_sub_temp + x[[jj_s]][ii_s, ]
    }
    means_sub[, ii_s] <- means_sub_temp / k
  }
  SSA <- rowSums(n * (means_gr - means_all)^2)
  SSS <- rowSums(k * (means_sub - means_all)^2)
  SST <- rowSums(sapply(lapply(x, function(x) (t(x) - means_all)^2), rowSums))
  SSE <- SST - SSA - SSS
  f_point <- (SSA / (k - 1)) / (SSE / ((n - 1) * (k - 1)))
  f_point <- f_point[is.finite(f_point)]
  f_point <- ifelse(f_point < .Machine$double.eps, 0, f_point)
  if (plot) {
    plot(f_point, type = type, ylab = ylab, main = main, ...)
  }
  if (values) {
    return(f_point)
  }
}

fpp_test_wiele_grup <- function(x, n_perm = 1000, n_boot = 1000,
                                parallel = FALSE, n_cores = NULL, multi_gen = FALSE) {
  if (parallel)
  {
    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
    l <- length(x)

    cn <- cn_value(x)
    dn_en_temp <- dn_en_value(x)
    dn <- dn_en_temp[1]
    en <- dn_en_temp[2]

    # test P1
    out <- foreach (i = 1:n_perm, .combine = rbind,
                    .packages = c("MASS", "rmfanova")) %dopar%
      {
        x_perm_1 <- vector("list", l)
        for (ii in 1:l) x_perm_1[[ii]] <- matrix(0, nrow = n, ncol = p)
        for (j in 1:n) {
          permutation <- sample(1:l)
          for (jj in 1:l) {
            x_perm_1[[jj]][j, ] <- x[[permutation[jj]]][j, ]
          }
        }
        dn_en_temp <- dn_en_value(x_perm_1)
        c(cn_value(x_perm_1), dn_en_temp[1], dn_en_temp[2])
      }
    pvalue_perm_1_cn <- mean(out[, 1] > cn)
    pvalue_perm_1_dn <- mean(out[, 2] > dn)
    pvalue_perm_1_en <- mean(out[, 3] > en)

    # test P2
    out2 <- foreach (i = 1:n_perm, .combine = rbind,
                     .packages = c("MASS", "rmfanova")) %dopar%
      {
        x_perm_2 <- vector("list", l)
        x_perm_mat <- x[[1]]
        for (ii in 1:l) {
          x_perm_2[[ii]] <- matrix(0, nrow = n, ncol = p)
          if (ii != 1) x_perm_mat <- rbind(x_perm_mat, x[[ii]])
        }
        for (jj in seq_len(l)) {
          permutation <- sample(seq_len(nrow(x_perm_mat)), n)
          x_perm_2[[jj]] <- x_perm_mat[permutation, ]
          x_perm_mat <- x_perm_mat[-permutation, ]
        }
        dn_en_temp <- dn_en_value(x_perm_2)
        c(cn_value(x_perm_2), dn_en_temp[1], dn_en_temp[2])
      }
    pvalue_perm_2_cn <- mean(out2[, 1] > cn)
    pvalue_perm_2_dn <- mean(out2[, 2] > dn)
    pvalue_perm_2_en <- mean(out2[, 3] > en)

    # test B1
    means_gr_raw <- sapply(x, colMeans)
    means_all_raw <- rowMeans(means_gr_raw)
    out3 <- foreach (i = 1:n_boot, .combine = rbind,
                     .packages = c("MASS", "rmfanova")) %dopar%
      {
        x_boot_1 <- vector("list", l)
        bootstrap <- floor(runif(n) * (n - 1)) + 1
        for (ii in 1:l) {
          x_boot_1[[ii]] <- x[[ii]][bootstrap, ]
        }
        dn_en_temp <- dn_en_value_boot(x_boot_1, means_gr_raw, means_all_raw)
        c(cn_value_boot(x_boot_1, means_gr_raw, means_all_raw), dn_en_temp[1], dn_en_temp[2])
      }
    pvalue_boot_1_cn <- mean(out3[, 1] > cn)
    pvalue_boot_1_dn <- mean(out3[, 2] > dn)
    pvalue_boot_1_en <- mean(out3[, 3] > en)

    # test B2
    means_gr_raw <- sapply(x, colMeans)
    out4 <- foreach (i = 1:n_boot, .combine = rbind,
                     .packages = c("MASS", "rmfanova")) %dopar%
      {
        x_boot_2 <- vector("list", l)
        for (ii in 1:l) {
          xii <- x[[ii]] - matrix(means_gr_raw[, ii], nrow = n, ncol = p, byrow = TRUE)
          x_boot_2[[ii]] <- xii[floor(runif(n) * (n - 1)) + 1, ]
        }
        dn_en_temp <- dn_en_value(x_boot_2)
        c(cn_value(x_boot_2), dn_en_temp[1], dn_en_temp[2])
      }
    pvalue_boot_2_cn <- mean(out4[, 1] > cn)
    pvalue_boot_2_dn <- mean(out4[, 2] > dn)
    pvalue_boot_2_en <- mean(out4[, 3] > en)

    # test B3
    x_matrix <- x[[1]]
    for (i_x_m in 2:l) {
      x_matrix <- cbind(x_matrix, x[[i_x_m]])
    }
    CC <- var(x_matrix)
    if (multi_gen == TRUE)
    {
      # proc_gauss <- MASS::mvrnorm(n = n * n_boot, rep(0, nrow(CC)), CC)
      out5 <- foreach (i_pb = 1:n_boot, .combine = rbind,
                       .packages = c("MASS", "rmfanova")) %dopar%
        {
          proc_gauss <- MASS::mvrnorm(n = n, rep(0, nrow(CC)), CC)
          # obserwacje_i_pb <- ((i_pb - 1) * n + 1):(i_pb * n)
          proc_gauss_list <- vector("list", l)
          for (i_p_g in seq_len(l)) {
            # proc_gauss_list[[i_p_g]] <- proc_gauss[obserwacje_i_pb, ((i_p_g - 1) * p + 1):(i_p_g * p)]
            proc_gauss_list[[i_p_g]] <- proc_gauss[, ((i_p_g - 1) * p + 1):(i_p_g * p)]
          }
          dn_en_temp <- dn_en_value(proc_gauss_list)
          c(cn_value(proc_gauss_list), dn_en_temp[1], dn_en_temp[2])
        }
    }
    else
    {
      proc_gauss <- MASS::mvrnorm(n = n * n_boot, rep(0, nrow(CC)), CC)
      out5 <- foreach (i_pb = seq_len(n_boot), .combine = rbind,
                       .packages = c("MASS", "rmfanova")) %dopar%
        {
          obserwacje_i_pb <- ((i_pb - 1) * n + 1):(i_pb * n)
          proc_gauss_list <- vector("list", l)
          for (i_p_g in seq_len(l)) {
            proc_gauss_list[[i_p_g]] <- proc_gauss[obserwacje_i_pb, ((i_p_g - 1) * p + 1):(i_p_g * p)]
          }
          dn_en_temp <- dn_en_value(proc_gauss_list)
          c(cn_value(proc_gauss_list), dn_en_temp[1], dn_en_temp[2])
        }
    }
    pvalue_pb_cn <- mean(out5[, 1] > cn)
    pvalue_pb_dn <- mean(out5[, 2] > dn)
    pvalue_pb_en <- mean(out5[, 3] > en)
  }
  else
  {
    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
    l <- length(x)

    cn <- cn_value(x)
    dn_en_temp <- dn_en_value(x)
    dn <- dn_en_temp[1]
    en <- dn_en_temp[2]

    # test P1
    cn_perm_1 <- numeric(n_perm)
    dn_perm_1 <- numeric(n_perm)
    en_perm_1 <- numeric(n_perm)
    for (i in 1:n_perm) {
      x_perm_1 <- vector("list", l)
      for (ii in 1:l) x_perm_1[[ii]] <- matrix(0, nrow = n, ncol = p)
      for (j in 1:n) {
        permutation <- sample(1:l)
        for (jj in 1:l) {
          x_perm_1[[jj]][j, ] <- x[[permutation[jj]]][j, ]
        }
      }
      cn_perm_1[i] <- cn_value(x_perm_1)
      dn_en_temp <- dn_en_value(x_perm_1)
      dn_perm_1[i] <- dn_en_temp[1]
      en_perm_1[i] <- dn_en_temp[2]
    }
    pvalue_perm_1_cn <- mean(cn_perm_1 > cn)
    pvalue_perm_1_dn <- mean(dn_perm_1 > dn)
    pvalue_perm_1_en <- mean(en_perm_1 > en)

    # test P2
    cn_perm_2 <- numeric(n_perm)
    dn_perm_2 <- numeric(n_perm)
    en_perm_2 <- numeric(n_perm)
    for (i in 1:n_perm) {
      x_perm_2 <- vector("list", l)
      x_perm_mat <- x[[1]]
      for (ii in 1:l) {
        x_perm_2[[ii]] <- matrix(0, nrow = n, ncol = p)
        if (ii != 1) x_perm_mat <- rbind(x_perm_mat, x[[ii]])
      }
      for (jj in seq_len(l)) {
        permutation <- sample(seq_len(nrow(x_perm_mat)), n)
        x_perm_2[[jj]] <- x_perm_mat[permutation, ]
        x_perm_mat <- x_perm_mat[-permutation, ]
      }
      cn_perm_2[i] <- cn_value(x_perm_2)
      dn_en_temp <- dn_en_value(x_perm_2)
      dn_perm_2[i] <- dn_en_temp[1]
      en_perm_2[i] <- dn_en_temp[2]
    }
    pvalue_perm_2_cn <- mean(cn_perm_2 > cn)
    pvalue_perm_2_dn <- mean(dn_perm_2 > dn)
    pvalue_perm_2_en <- mean(en_perm_2 > en)

    # test B1
    means_gr_raw <- sapply(x, colMeans)
    means_all_raw <- rowMeans(means_gr_raw)
    cn_boot_1 <- numeric(n_boot)
    dn_boot_1 <- numeric(n_boot)
    en_boot_1 <- numeric(n_boot)
    for (i in 1:n_boot) {
      x_boot_1 <- vector("list", l)
      bootstrap <- floor(runif(n) * (n - 1)) + 1
      for (ii in 1:l) {
        x_boot_1[[ii]] <- x[[ii]][bootstrap, ]
      }
      cn_boot_1[i] <- cn_value_boot(x_boot_1, means_gr_raw, means_all_raw)
      dn_en_temp <- dn_en_value_boot(x_boot_1, means_gr_raw, means_all_raw)
      dn_boot_1[i] <- dn_en_temp[1]
      en_boot_1[i] <- dn_en_temp[2]
    }
    pvalue_boot_1_cn <- mean(cn_boot_1 > cn)
    pvalue_boot_1_dn <- mean(dn_boot_1 > dn)
    pvalue_boot_1_en <- mean(en_boot_1 > en)

    # test B2
    means_gr_raw <- sapply(x, colMeans)
    cn_boot_2 <- numeric(n_boot)
    dn_boot_2 <- numeric(n_boot)
    en_boot_2 <- numeric(n_boot)
    for (i in 1:n_boot) {
      x_boot_2 <- vector("list", l)
      for (ii in 1:l) {
        xii <- x[[ii]] - matrix(means_gr_raw[, ii], nrow = n, ncol = p, byrow = TRUE)
        x_boot_2[[ii]] <- xii[floor(runif(n) * (n - 1)) + 1, ]
      }
      cn_boot_2[i] <- cn_value(x_boot_2)
      dn_en_temp <- dn_en_value(x_boot_2)
      dn_boot_2[i] <- dn_en_temp[1]
      en_boot_2[i] <- dn_en_temp[2]
    }
    pvalue_boot_2_cn <- mean(cn_boot_2 > cn)
    pvalue_boot_2_dn <- mean(dn_boot_2 > dn)
    pvalue_boot_2_en <- mean(en_boot_2 > en)

    # test B3
    x_matrix <- x[[1]]
    for (i_x_m in 2:l) {
      x_matrix <- cbind(x_matrix, x[[i_x_m]])
    }
    CC <- var(x_matrix)
    cn_pb <- numeric(n_boot)
    dn_pb <- numeric(n_boot)
    en_pb <- numeric(n_boot)
    if (multi_gen == TRUE)
    {
      # proc_gauss <- MASS::mvrnorm(n = n * n_boot, rep(0, nrow(CC)), CC)
      for (i_pb in seq_len(n_boot)) {
        proc_gauss <- MASS::mvrnorm(n = n, rep(0, nrow(CC)), CC)
        # obserwacje_i_pb <- ((i_pb - 1) * n + 1):(i_pb * n)
        proc_gauss_list <- vector("list", l)
        for (i_p_g in seq_len(l)) {
          # proc_gauss_list[[i_p_g]] <- proc_gauss[obserwacje_i_pb, ((i_p_g - 1) * p + 1):(i_p_g * p)]
          proc_gauss_list[[i_p_g]] <- proc_gauss[, ((i_p_g - 1) * p + 1):(i_p_g * p)]
        }
        cn_pb[i_pb] <- cn_value(proc_gauss_list)
        dn_en_temp <- dn_en_value(proc_gauss_list)
        dn_pb[i_pb] <- dn_en_temp[1]
        en_pb[i_pb] <- dn_en_temp[2]
      }
    }
    else
    {
      proc_gauss <- MASS::mvrnorm(n = n * n_boot, rep(0, nrow(CC)), CC)
      for (i_pb in seq_len(n_boot)) {
        obserwacje_i_pb <- ((i_pb - 1) * n + 1):(i_pb * n)
        proc_gauss_list <- vector("list", l)
        for (i_p_g in seq_len(l)) {
          proc_gauss_list[[i_p_g]] <- proc_gauss[obserwacje_i_pb, ((i_p_g - 1) * p + 1):(i_p_g * p)]
        }
        cn_pb[i_pb] <- cn_value(proc_gauss_list)
        dn_en_temp <- dn_en_value(proc_gauss_list)
        dn_pb[i_pb] <- dn_en_temp[1]
        en_pb[i_pb] <- dn_en_temp[2]
      }
    }
    pvalue_pb_cn <- mean(cn_pb > cn)
    pvalue_pb_dn <- mean(dn_pb > dn)
    pvalue_pb_en <- mean(en_pb > en)
  }
  return(c(pvalue_perm_1_cn, pvalue_perm_2_cn, pvalue_boot_1_cn, pvalue_boot_2_cn, pvalue_pb_cn,
           pvalue_perm_1_dn, pvalue_perm_2_dn, pvalue_boot_1_dn, pvalue_boot_2_dn, pvalue_pb_dn,
           pvalue_perm_1_en, pvalue_perm_2_en, pvalue_boot_1_en, pvalue_boot_2_en, pvalue_pb_en))
}

#' Repeated measures functional analysis of variance
#'
#' The function \code{rmfanova()} calculates the tests based on three test statistics
#' \eqn{\mathcal{C}_n}, \eqn{\mathcal{D}_n}, and \eqn{\mathcal{E}_n} for the problem of
#' comparing \eqn{\ell}-samples of repeated measures for functional data. The tests are based on
#' five resampling methods, i.e., two permutation and three bootstrap ones. The overall and local
#' hypotheses are considered.
#'
#' @param x a list of length \eqn{\ell} with elements being \eqn{n\times p} matrices of data
#' corresponding to \eqn{n} functional observations measured in \eqn{p} design time points under given
#' experimental conditions.
#' @param method the correction method to be used for pairwise comparisons. Options are \code{"bonferroni"}
#' (default) and those given in the vector \code{p.adjust.methods} (as for the \code{p.adjust()} function).
#' @param n_perm a number of permutation replicates. The default is 1000.
#' @param n_boot a number of bootstrap replicates. The default is 1000.
#' @param parallel a logical indicating of whether to use parallel computing. The default is \code{FALSE.}
#' @param n_cores if \code{parallel = TRUE}, a number of processes used in parallel computation.
#' Its default value (\code{NULL}) means that it will be equal to a number of cores of a computer used.
#' @param multi_gen a logical indicating of whether to use separate multiple generations of Gaussian processes
#' for the parametric bootstrap tests. The default is FALSE, which means that the processes will be
#' generated once in a big matrix. This method is much faster, but for larger \eqn{n} and \eqn{p}
#' the generated data can be too large for RAM. In such a case, we suggest using separate generation
#' (\code{multi_gen = TRUE}), which is slower, but possible to calculate.
#'
#' @details The function \code{rmfanova()} concerns the tests for the functional repeated measures analysis problem.
#' The details are presented in Kurylo and Smaga (2023), where in particular, some recommendations for using tests are given.
#' Here we present only some summary of the problem and its solutions implemented in the package.
#'
#' We have \eqn{n} subjects subjected to \eqn{\ell\geq 2} (possibly) different conditions.
#' The results of the experiments are functional observations. Let the subjects be represented
#' by a functional sample consisting of independent stochastic processes \eqn{Y_1,\dots,Y_n} defined on the
#' interval \eqn{[0,\ell]}, which satisfy the following model proposed by Martinez-Camblor and Corral (2011):
#' \deqn{Y_j(t)=\mu(t)+e_j(t),\ j=1,\dots,n,\ t\in[0,\ell],}
#' where \eqn{\mu} is a fixed mean function, and \eqn{e_j} is a random process with zero mean function.
#' In this notation, \eqn{t\in[0,1]} corresponds to the first experimental condition, \eqn{t\in[1,2]}
#' to the second, and so on. Thus, in this model, we ignore the possible time periods between repetitions
#' of the experiment, but this does not mean that they do not exist. We are interested in testing the equality
#' of \eqn{\ell} mean functions corresponding to experimental conditions; namely, the global null hypothesis is as follows:
#' \deqn{\mathcal{H}_0:\mu(t)=\mu(t+1)=\dots=\mu(t+(\ell-1))\ \ \forall t\in[0,1].}
#' For the global null hypothesis \eqn{\mathcal{H}_0}, the tests given by Martinez-Camblor and Corral (2011)
#' used the pointwise sum of squares due to the hypothesis:
#' \deqn{\mathrm{SSA}_{point}(t)=n\sum_{i=1}^\ell(\bar{Y}_{i\cdot}(t)-\bar{Y}(t))^2,\ t\in[0,1],}
#' where \deqn{\bar{Y}_{i\cdot}(t)=n^{-1}\sum_{j=1}^nY_j(t+(i-1)),\ \bar{Y}(t)=N^{-1}\sum_{i=1}^\ell\sum_{j=1}^nY_j(t+(i-1)),}
#' \eqn{i=1,\dots,\ell}. In the package, it is calculated and drawn by the \code{pointwise_ssa_test_statistic()} function.
#' The other option is the following pointwise F-type test statistic proposed in Kurylo and Smaga (2023):
#' \deqn{F_{point}(t)=\frac{\mathrm{SSA}_{point}(t)/(\ell-1)}{\mathrm{SSR}_{point}(t)/((\ell-1)(n-1))},\ t\in[0,1],}
#' where \deqn{\mathrm{SSR}_{point}(t)=\sum_{i=1}^\ell\sum_{j=1}^n(Y_j(t+(i-1))-\bar{Y}_{i\cdot}(t)-\bar{Y}_{\cdot j}(t)+\bar{Y}(t))^2}
#' is the pointwise sum of squares due to residuals, and \deqn{\bar{Y}_{\cdot j}(t)=\ell^{-1}\sum_{i=1}^\ell Y_j(t+(i-1)),\ j=1,\dots,n.}
#' \eqn{F_{point}} is calculated and drawn by the \code{pointwise_f_test_statistic()} function.
#'
#' To obtain global test statistics for \eqn{\mathcal{H}_0}, Martinez-Camblor and Corral (2011) proposed the
#' following test statistic: \deqn{\mathcal{C}_n(\ell)=\int_0^1\mathrm{SSA}_{point}(t)dt.} On the other hand,
#' Kurylo and Smaga (2023) proposed the following two test statistics:
#' \deqn{\mathcal{D}_n(\ell)=\int_0^1F_{point}(t)dt,\quad\mathcal{E}_n(\ell)=\sup\limits_{t\in[0,1]}F_{point}(t).}
#' To construct the tests, five resampling strategies are proposed by Kurylo and Smaga (2023). For details, we refer
#' to this paper. Here we just note the two permutation tests and three bootstrap tests are denoted by P1, P2, B1, B2,
#' and B3 in the output of the \code{summary.rmfanova()} function.
#'
#' When \eqn{\ell>2}, by rejecting the global null hypothesis \eqn{\mathcal{H}_0}, we determine the presence of significant differences
#' in the mean functions corresponding to the experimental conditions. However, we do not know which conditions are
#' significantly different and which are not. To solve this problem, one needs to perform a post hoc analysis.
#' More precisely, we would like to test the family of hypotheses:
#' \deqn{\left\{\begin{array}{l}
#' \mathcal{H}_0^{rs}:\mu(t+(r-1))=\mu(t+(s-1))\ \forall t\in[0,1],\\
#' \mathcal{H}_1^{rs}:\mu(t+(r-1))\neq\mu(t+(s-1))\ \text{for some}\ t\in[0,1],\\
#' \end{array}\right.}
#' for \eqn{r,s=1,\dots,\ell}, \eqn{r\neq s}. These hypotheses are also named pairwise comparisons.
#' To test this family of local hypotheses, we propose the following procedure:
#'
#' 1. Test each of the hypotheses \eqn{\mathcal{H}_0^{rs}}  using the data for the \eqn{r}-th and \eqn{s}-th objects,
#' i.e., \eqn{Y_1(t),\dots,Y_n(t)} for \eqn{t\in[r-1,r]} and \eqn{t\in[s-1,s]} respectively, and the chosen test
#' from those presented above. Let \eqn{p_{rs}} denote the \eqn{p}-values obtained.
#'
#' 2. Make a final decision using the Bonferroni method, i.e., reject \eqn{\mathcal{H}_0^{rs}} if
#' \eqn{p_{rs}^{Bonf}\leq \alpha}, where \eqn{p_{rs}^{Bonf}=m\cdot p_{rs}} are the corrected \eqn{p}-values,
#' \eqn{\alpha} is the significance level and \eqn{m} is the number of null hypotheses considered.
#'
#' In the paper Kurylo and Smaga (2023), the Bonferroni method was used only. However, in the package,
#' there is a possibility to use other correction methods, which are available in the vector \code{p.adjust.methods}.
#'
#' The results of testing the global and local hypotheses are given separately in the output of the
#' \code{summary.rmfanova()} function for the convenience of the user.
#'
#' @return A list of class \code{rmfanova} containing the following 7 components:
#' \item{n}{a number \eqn{n} of functional observations.}
#' \item{p}{a number \eqn{p} of design time points.}
#' \item{l}{a number \eqn{\ell} of repeated samples.}
#' \item{method}{an argument \code{method}.}
#' \item{test_stat}{values of the test statistics \eqn{\mathcal{C}_n}, \eqn{\mathcal{D}_n}, and \eqn{\mathcal{E}_n}.}
#' \item{p_values}{p-values for the global null hypothesis.}
#' \item{p_values_pc}{p-values of the pairwise comparisons.}
#'
#' @references Martinez-Camblor P., Corral N. (2011) Repeated Measures Analysis for Functional Data.
#' Computational Statistics & Data Analysis 55, 3244–3256.
#'
#' Kurylo K., Smaga L. (2023) Functional repeated measures analysis of variance and its application. Preprint https://arxiv.org/abs/2306.03883
#'
#' Ramsay J.O., Silverman B.W. (2005) Functional Data Analysis, 2nd Edition. New York: Springer.
#'
#' Zhang J.T. (2013) Analysis of Variance for Functional Data. London: Chapman & Hall.
#'
#' @examples
#' # Some of the examples may run some time.
#' # preparation of the DTI data set, for details see Kurylo and Smaga (2023)
#' library(refund)
#' data(DTI)
#' # MS patients
#' DTI_ms <- DTI[DTI$case == 1, ]
#' miss_data <- c()
#' for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
#' DTI_ms <- DTI_ms[-miss_data, ]
#' DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
#' xx <- vector("list", 4)
#' for (i in 1:4) {
#'   xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
#' }
#' xx[[1]] <- xx[[1]][-14, ]
#' xx[[3]] <- xx[[3]][-14, ]
#' yy <- xx
#' for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
#' # data trajectories for four visits
#' oldpar <- par(mfrow = c(1, 4), mar = c(4, 4, 4, 0.1))
#' matplot(t(yy[[1]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
#'         main = "Visit 1", xaxt = "n", ylim = c(0.29, 0.73))
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' matplot(t(yy[[2]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
#'         main = "Visit 2", xaxt = "n", ylim = c(0.29, 0.73))
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' matplot(t(yy[[3]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
#'         main = "Visit 3", xaxt = "n", ylim = c(0.29, 0.73))
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' matplot(t(yy[[4]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
#'         main = "Visit 4", xaxt = "n", ylim = c(0.29, 0.73))
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' par(oldpar)
#' # sample mean functions
#' oldpar <- par(mfrow = c(1, 1), mar = c(4, 4, 2, 0.1))
#' pointwise_sample_mean_fun(yy, values = FALSE,
#'                           col = 1:4, xlab = "t", ylab = "FA", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' legend(x = 36, y = 0.64, legend = 1:4, lty = 1, col = 1:4, title = "Visit")
#' par(oldpar)
#' # pointwise SSA and F-type test statistics
#' oldpar <- par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1))
#' pointwise_ssa_test_statistic(yy, xlab = "t", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' pointwise_f_test_statistic(yy, xlab = "t", xaxt = "n")
#' axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
#' par(oldpar)
#' \donttest{
#' # testing without parallel computing and multiple generation of Gaussian processes
#' res <- rmfanova(yy)
#' summary(res, digits = 3)
#' # testing without parallel computing and with multiple generation of Gaussian processes
#' res <- rmfanova(yy, multi_gen = TRUE)
#' summary(res, digits = 3)
#' # testing with parallel computing and without multiple generation of Gaussian processes
#' res <- rmfanova(yy, parallel = TRUE, n_cores = 2)
#' summary(res, digits = 3)
#' # testing with parallel computing and with multiple generation of Gaussian processes
#' res <- rmfanova(yy, parallel = TRUE, multi_gen = TRUE, n_cores = 2)
#' summary(res, digits = 3)}
#' \dontshow{
#' yy_2 <- list(yy[[1]][1:5, 1:2],
#'              yy[[2]][1:5, 1:2],
#'              yy[[3]][1:5, 1:2])
#' # testing without parallel computing and multiple generation of Gaussian processes
#' res <- rmfanova(yy_2, n_perm = 10, n_boot = 10)
#' summary(res, digits = 3)
#' # testing without parallel computing and with multiple generation of Gaussian processes
#' res <- rmfanova(yy_2, n_perm = 10, n_boot = 10, multi_gen = TRUE)
#' summary(res, digits = 3)
#' # testing with parallel computing and without multiple generation of Gaussian processes
#' res <- rmfanova(yy_2, n_perm = 2, n_boot = 2, parallel = TRUE, n_cores = 2)
#' summary(res, digits = 3)}
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import refund
#' @import MASS
#' @importFrom graphics matplot
#' @importFrom stats p.adjust runif var
#'
#' @export
rmfanova <- function(x, method = "bonferroni",
                     n_perm = 1000, n_boot = 1000,
                     parallel = FALSE, n_cores = NULL, multi_gen = FALSE) {
  if (is.list(x) == FALSE) {
    stop("x should be a list")
  }
  if (is.numeric(x[[1]]) == FALSE) {
    stop("x should have numeric values")
  }
  if (n_boot %% 1 != 0) {
    stop("n_boot should be integer")
  }
  if (n_perm %% 1 != 0) {
    stop("n_perm should be integer")
  }
  if (parallel)
  {
    # #' @importFrom utils installed.packages
    # if (!("doParallel" %in% rownames(installed.packages()))) {
    #   stop("Please install package 'doParallel'")
    # }
    # require(foreach, quietly = TRUE)
    requireNamespace("doParallel", quietly = TRUE)
    requireNamespace("foreach", quietly = TRUE)
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores()
    } else if (!(n_cores > 1)) {
      stop("n_cores should be greater than 1")
    }
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)

    ll <- length(x)
    # main tests
    p_val <- fpp_test_wiele_grup(x, n_perm = n_perm, n_boot = n_boot,
                                 parallel = parallel, n_cores = n_cores,
                                 multi_gen = multi_gen)
    p_val <- t(as.data.frame(p_val))
    rownames(p_val) <- "1"
    colnames(p_val) <- c("Cn_P1", "Cn_P2", "Cn_B1", "Cn_B2", "Cn_B3",
                         "Dn_P1", "Dn_P2", "Dn_B1", "Dn_B2", "Dn_B3",
                         "En_P1", "En_P2", "En_B1", "En_B2", "En_B3")
    # pairwise comparisons
    if (ll == 2) {
      p_val_pc <- NULL
    } else if (ll > 2) {
      p_val_pc <- matrix(NA, nrow = choose(ll, 2), ncol = 15)
      i_temp <- 1
      for (i in 1:(ll - 1)) {
        for (j in (i + 1):ll) {
          p_val_pc[i_temp, ] <- fpp_test_wiele_grup(x[c(i, j)], n_perm = n_perm, n_boot = n_boot,
                                                    parallel = parallel, n_cores = n_cores,
                                                    multi_gen = multi_gen)
          i_temp <- i_temp + 1
        }
      }
      p_val_pc <- apply(p_val_pc, 2, p.adjust, method = method)
      p_val_pc <- as.data.frame(p_val_pc)
      colnames(p_val_pc) <- c("Cn_P1", "Cn_P2", "Cn_B1", "Cn_B2", "Cn_B3",
                              "Dn_P1", "Dn_P2", "Dn_B1", "Dn_B2", "Dn_B3",
                              "En_P1", "En_P2", "En_B1", "En_B2", "En_B3")
      i_temp <- 1
      for (i in 1:(ll - 1)) {
        for (j in (i + 1):ll) {
          rownames(p_val_pc)[i_temp] <- paste(i, "-", j)
          i_temp <- i_temp + 1
        }
      }
    }
    res <- list(n = nrow(x[[1]]), p = ncol(x[[1]]), l = ll, method = method,
                test_stat = data.frame(Cn = cn_value(x),
                                       Dn = dn_en_value(x)[1],
                                       En = dn_en_value(x)[2]),
                p_values = p_val,
                p_values_pc = p_val_pc)
    class(res) <- "rmfanova"
  }
  else
  {
    ll <- length(x)
    # main tests
    p_val <- fpp_test_wiele_grup(x, n_perm = n_perm, n_boot = n_boot,
                                 parallel = parallel, n_cores = n_cores,
                                 multi_gen = multi_gen)
    p_val <- t(as.data.frame(p_val))
    rownames(p_val) <- "1"
    colnames(p_val) <- c("Cn_P1", "Cn_P2", "Cn_B1", "Cn_B2", "Cn_B3",
                         "Dn_P1", "Dn_P2", "Dn_B1", "Dn_B2", "Dn_B3",
                         "En_P1", "En_P2", "En_B1", "En_B2", "En_B3")
    # pairwise comparisons
    if (ll == 2) {
      p_val_pc <- NULL
    } else if (ll > 2) {
      p_val_pc <- matrix(NA, nrow = choose(ll, 2), ncol = 15)
      i_temp <- 1
      for (i in 1:(ll - 1)) {
        for (j in (i + 1):ll) {
          p_val_pc[i_temp, ] <- fpp_test_wiele_grup(x[c(i, j)], n_perm = n_perm, n_boot = n_boot,
                                                    parallel = parallel, n_cores = n_cores,
                                                    multi_gen = multi_gen)
          i_temp <- i_temp + 1
        }
      }
      p_val_pc <- apply(p_val_pc, 2, p.adjust, method = method)
      p_val_pc <- as.data.frame(p_val_pc)
      colnames(p_val_pc) <- c("Cn_P1", "Cn_P2", "Cn_B1", "Cn_B2", "Cn_B3",
                              "Dn_P1", "Dn_P2", "Dn_B1", "Dn_B2", "Dn_B3",
                              "En_P1", "En_P2", "En_B1", "En_B2", "En_B3")
      i_temp <- 1
      for (i in 1:(ll - 1)) {
        for (j in (i + 1):ll) {
          rownames(p_val_pc)[i_temp] <- paste(i, "-", j)
          i_temp <- i_temp + 1
        }
      }
    }
    res <- list(n = nrow(x[[1]]), p = ncol(x[[1]]), l = ll, method = method,
                test_stat = data.frame(Cn = cn_value(x),
                                       Dn = dn_en_value(x)[1],
                                       En = dn_en_value(x)[2]),
                p_values = p_val,
                p_values_pc = p_val_pc)
    class(res) <- "rmfanova"
  }
  return(res)
}

#' Print "rmfanova" object
#'
#' Prints the summary of the repeated measures functional analysis of variance.
#'
#' @param object a "rmfanova" object.
#' @param ... integer indicating the number of decimal places to be used to present the numerical results.
#' It can be named \code{digits} as in the \code{round()} function (see examples).
#'
#' @details The function prints out the information about the number of samples \eqn{\ell},
#' number of observations \eqn{n}, number of design time points \eqn{p},
#' adjustment method for pairwise comparison tests (if \eqn{\ell>2}), test statistics,
#' and p-values of tests performed by the \code{rmfanova()} function.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' # Some of the examples may run some time.
#' # preparation of the DTI data set, for details see Kurylo and Smaga (2023)
#' library(refund)
#' data(DTI)
#' # MS patients
#' DTI_ms <- DTI[DTI$case == 1, ]
#' miss_data <- c()
#' for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
#' DTI_ms <- DTI_ms[-miss_data, ]
#' DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
#' xx <- vector("list", 4)
#' for (i in 1:4) {
#'   xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
#' }
#' xx[[1]] <- xx[[1]][-14, ]
#' xx[[3]] <- xx[[3]][-14, ]
#' yy <- xx
#' for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
#' \donttest{
#' # testing without parallel computing and multiple generation of Gaussian processes
#' res <- rmfanova(yy)
#' summary(res, digits = 3)}
#'
#' @export
summary.rmfanova <- function(object, ...) {
  cat("#--- Repeated measures funtional analysis of variance -----------------------------#", "\n", "\n")
  cat("- Number of samples:", object$l, "\n")
  cat("- Number of observations:", object$n, "\n")
  cat("- Number of design time points:", object$p, "\n")
  if (object$l > 2) {
    cat("- Adjustment method for pairwise comparison tests:", object$method)
  }

  cat("\n", "\n")
  cat("#--- Overall test statistics ------------------------------------------------------#", "\n")
  print(round(object$test_stat, ...))
  cat("\n")
  cat("#--- Overall p-values -------------------------------------------------------------#", "\n")
  print(round(object$p_values, ...))
  if (object$l > 2) {
    cat("\n")
    cat("#--- Pairwise comparison p-values -------------------------------------------------#", "\n")
    print(round(object$p_values_pc, ...))
  }
  cat("#----------------------------------------------------------------------------------#", "\n")
}
