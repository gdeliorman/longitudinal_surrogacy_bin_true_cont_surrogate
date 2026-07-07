rm(list = ls())

## libraries
suppressPackageStartupMessages({
  library(MCMCglmm)
  library(lme4)
  library(readr)
  library(dplyr)
  library(tidyr)
})


## data
my_data   <- read.csv("/Users/gokcedeliorman/Downloads/data_binary_cgi_panss.csv")
wide_data <- read.csv("/Users/gokcedeliorman/Downloads/wide_data_binary_cgi.csv")

colnames(wide_data)
str(wide_data)
head(wide_data)
colSums(is.na(wide_data))

df_joint <- wide_data %>%
  mutate(
    id    = as.factor(id),
    TREAT = as.factor(TREAT),
    time  = as.factor(time),
    XTIME_cont = as.numeric(XTIME_cont),
    CGI   = factor(as.numeric(response_log_CGI), levels = c(0, 1)),
    PANSS = as.numeric(response_log_PANSS),
    base_CGI_log   = as.numeric(baseline_value_CGI),
    base_PANSS_log = as.numeric(baseline_value_PANSS),
    base_CGI_bin = as.factor(base_CGI_bin),
    binary_CGI =as.factor(binary_CGI),
    base_CGI_bin = as.numeric(as.character(base_CGI_bin))
  ) 

str(df_joint)
df_joint$TREAT <- factor(df_joint$TREAT, levels = c("0", "1"))


df_joint$id_t0 <- factor(ifelse(df_joint$TREAT == "0", as.character(df_joint$id), NA))
df_joint$id_t1 <- factor(ifelse(df_joint$TREAT == "1", as.character(df_joint$id), NA))


#trait = 1 → PANSS (S)
#trait = 2 → CGI (T)

fam <- c("gaussian", "threshold")

f_random_ri <- ~ us(trait):id_t0 + us(trait):id_t1
f_random_rs <- ~ us(trait + trait:XTIME_cont):id_t0 +
  us(trait + trait:XTIME_cont):id_t1

f_shared_random_ri <- ~ id_t0 + id_t1
f_shared_random_rs <- ~ us(1 + XTIME_cont):id_t0 +
  us(1 + XTIME_cont):id_t1

f_fixed <- cbind(PANSS, CGI) ~ 0 + trait:time:TREAT +
  at.level(trait, 1):base_PANSS_log:TREAT +
  at.level(trait, 2):base_CGI_bin:TREAT


# f_fixed_1 <- cbind(PANSS, CGI) ~ 0 + trait:time+ trait:time:TREAT +
#   at.level(trait, 1):base_PANSS_log +
#   at.level(trait, 2):base_CGI_bin
# 
# 
# f_fixed_2 <- cbind(PANSS, CGI) ~ 0 + trait +
#    at.level(trait, 1):(base_PANSS_log + time*TREAT) +
#    at.level(trait, 2):(base_CGI_bin + time*TREAT)

# same as fixed_2 : f_fixed_3 <- cbind(PANSS, CGI) ~ trait - 1 +
#   at.level(trait, 1):(base_PANSS_log + time*TREAT) +
#   at.level(trait, 2):(base_CGI_bin + time*TREAT)


## weak proper priors
rcov_trt <- ~
  idh(trait:at.level(TREAT, "0")):units +
  idh(trait:at.level(TREAT, "1")):units

R_prior_2blocks <- list(
  R1 = list(V = diag(c(0.02, 1)), nu = 2, fix = 2),  # TREAT=0
  R2 = list(V = diag(c(0.02, 1)), nu = 2, fix = 2)   # TREAT=1
)

prior_ri <- list(
  G = list(
    G1 = list(V = diag(2), nu = 2.002),
    G2 = list(V = diag(2), nu = 2.002)
  ),
  R = R_prior_2blocks
)

prior_rs <- list(
  G = list(
    G1 = list(V = diag(4), nu = 4.002),
    G2 = list(V = diag(4), nu = 4.002)
  ),
  R = R_prior_2blocks
)

## MCMC settings to test
nitt_use   <- 80
burnin_use <- 20
thin_use   <- 3

## MCMC real settings
nitt_use   <- 800000
burnin_use <- 200000
thin_use   <- 300

start_time <- Sys.time()
set.seed(2002)

m_ri <- MCMCglmm(
  fixed   = f_fixed,
  random  = f_random_ri,
  rcov    = rcov_trt,
  family  = fam,
  data    = df_joint,
  prior   = prior_ri,
  nitt    = nitt_use,
  burnin  = burnin_use,
  thin    = thin_use,
  verbose = TRUE
)

m_rs <- MCMCglmm(
  fixed   = f_fixed,
  random  = f_random_rs,
  rcov    = rcov_trt,
  family  = fam,
  data    = df_joint,
  prior   = prior_rs,
  nitt    = nitt_use,
  burnin  = burnin_use,
  thin    = thin_use,
  verbose = TRUE
)


m_ri_shared <- MCMCglmm(
  fixed   = f_fixed,
  random  = f_shared_random_ri,
  rcov    = rcov_trt,
  family  = fam,
  data    = df_joint,
  prior   = prior_ri,
  nitt    = nitt_use,
  burnin  = burnin_use,
  thin    = thin_use,
  verbose = TRUE
)

m_rs_shared <- MCMCglmm(
  fixed   = f_fixed,
  random  = f_shared_random_rs,
  rcov    = rcov_trt,
  family  = fam,
  data    = df_joint,
  prior   = prior_rs,
  nitt    = nitt_use,
  burnin  = burnin_use,
  thin    = thin_use,
  verbose = TRUE
)


end_time <- Sys.time()



## ------------------------------------------------------------
## helper to extract D matrix
## ------------------------------------------------------------
get_R_4x4 <- function(fit) {
  m <- colMeans(fit$VCV)
  
  R0 <- matrix(
    c(
      m['traitCGI:at.level(TREAT, "0").units'], 0,
      0, m['traitPANSS:at.level(TREAT, "0").units']
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("T0", "S0"), c("T0", "S0"))
  )
  
  R1 <- matrix(
    c(
      m['traitCGI:at.level(TREAT, "1").units'], 0,
      0, m['traitPANSS:at.level(TREAT, "1").units']
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("T1", "S1"), c("T1", "S1"))
  )
  
  R <- matrix(
    0, 4, 4,
    dimnames = list(c("T0", "T1", "S0", "S1"),
                    c("T0", "T1", "S0", "S1"))
  )
  
  R["T0", "T0"] <- m['traitCGI:at.level(TREAT, "0").units']
  R["T1", "T1"] <- m['traitCGI:at.level(TREAT, "1").units']
  R["S0", "S0"] <- m['traitPANSS:at.level(TREAT, "0").units']
  R["S1", "S1"] <- m['traitPANSS:at.level(TREAT, "1").units']
  
  list(R0 = R0, R1 = R1, R = R)
}


get_R_4x4(m_ri)
get_R_4x4(m_rs)

R_ri<-get_R_4x4(m_ri)$R
R_rs<-get_R_4x4(m_rs)$R


get_D <- function(fit,
                  blocks,
                  labels_by_block,
                  names_by_block = labels_by_block,
                  display_order = NULL) {
  m <- colMeans(fit$VCV)
  
  get_block <- function(block, labels, pretty_names) {
    D <- matrix(
      NA_real_,
      nrow = length(labels),
      ncol = length(labels),
      dimnames = list(pretty_names, pretty_names)
    )
    
    for (i in seq_along(labels)) {
      for (j in seq_along(labels)) {
        key_ij <- paste0(labels[i], ":", labels[j], ".", block)
        key_ji <- paste0(labels[j], ":", labels[i], ".", block)
        
        if (key_ij %in% names(m)) {
          D[i, j] <- m[[key_ij]]
        } else if (key_ji %in% names(m)) {
          D[i, j] <- m[[key_ji]]
        } else {
          stop("Missing key in VCV: ", key_ij)
        }
      }
    }
    
    D
  }
  
  block_list <- Map(get_block, blocks, labels_by_block, names_by_block)
  names(block_list) <- blocks
  
  total_dim <- sum(vapply(block_list, nrow, integer(1)))
  all_names <- unlist(lapply(block_list, rownames), use.names = FALSE)
  
  D_full <- matrix(
    0,
    nrow = total_dim,
    ncol = total_dim,
    dimnames = list(all_names, all_names)
  )
  
  idx <- 1
  for (B in block_list) {
    k <- nrow(B)
    D_full[idx:(idx + k - 1), idx:(idx + k - 1)] <- B
    idx <- idx + k
  }
  
  if (!is.null(display_order)) {
    D_full <- D_full[display_order, display_order, drop = FALSE]
  }
  
  list(
    blocks = block_list,
    D = D_full
  )
}

D_ri_out <- get_D(
  m_ri,
  blocks = c("id_t0", "id_t1"),
  labels_by_block = list(
    c("traitPANSS", "traitCGI"),
    c("traitPANSS", "traitCGI")
  ),
  names_by_block = list(
    c("S0", "T0"),
    c("S1", "T1")
  ),
  display_order = c("T0", "T1", "S0", "S1")
)


D_ri_out

D_rs_out <- get_D(
  m_rs,
  blocks = c("id_t0", "id_t1"),
  labels_by_block = list(
    c("traitPANSS", "traitCGI", "traitPANSS:XTIME_cont", "traitCGI:XTIME_cont"),
    c("traitPANSS", "traitCGI", "traitPANSS:XTIME_cont", "traitCGI:XTIME_cont")
  ),
  names_by_block = list(
    c("S0_int", "T0_int", "S0_slope", "T0_slope"),
    c("S1_int", "T1_int", "S1_slope", "T1_slope")
  ),
  display_order = c("T0_int", "T0_slope", "T1_int", "T1_slope",
                    "S0_int", "S0_slope", "S1_int", "S1_slope")
)


D_rs_out

D_ri<- D_ri_out$D
D_rs<- D_rs_out$D

summarise_posterior_matrix <- function(draws) {
  out <- t(apply(draws, 2, function(x) {
    c(
      PostMean = mean(x),
      CrI_Lower = unname(quantile(x, 0.025)),
      CrI_Upper = unname(quantile(x, 0.975))
    )
  }))
  
  data.frame(
    Parameter = rownames(out),
    PostMean = out[, "PostMean"],
    CrI_Lower = out[, "CrI_Lower"],
    CrI_Upper = out[, "CrI_Upper"],
    row.names = NULL,
    check.names = FALSE
  )
}

fixed_summary_ri <- summarise_posterior_matrix(m_ri$Sol)
fixed_summary_rs <- summarise_posterior_matrix(m_rs$Sol)
#fixed_summary_ri2 <- summarise_posterior_matrix(m_ri2$Sol)
#fixed_summary_rs2 <- summarise_posterior_matrix(m_rs2$Sol)

vcv_summary_ri <- summarise_posterior_matrix(m_ri$VCV)
vcv_summary_rs <- summarise_posterior_matrix(m_rs$VCV)

model_fit_summary <- data.frame(
  Model = c("m_ri", "m_rs"),
  PosteriorSamples_Sol = c(nrow(m_ri$Sol), nrow(m_rs$Sol)),
  PosteriorSamples_VCV = c(nrow(m_ri$VCV), nrow(m_rs$VCV)),
  row.names = NULL
)



parameter_summary <- rbind(
  cbind(Model = "m_ri", Group = "Fixed", fixed_summary_ri, row.names = NULL),
  cbind(Model = "m_ri", Group = "VCV",   vcv_summary_ri, row.names = NULL),
  cbind(Model = "m_rs", Group = "Fixed", fixed_summary_rs, row.names = NULL),
  cbind(Model = "m_rs", Group = "VCV",   vcv_summary_rs, row.names = NULL)
)


## ------------------------------------------------------------
## save outputs
## ------------------------------------------------------------

saveRDS(m_ri, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/m_ri.rds")
saveRDS(m_rs, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/m_rs.rds")

saveRDS(m_rs_shared, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/m_rs_shared.rds")
saveRDS(m_ri_shared, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/m_ri_shared.rds")

write.csv(D_ri, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/D_ri.csv", row.names = TRUE)
write.csv(D_rs, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/D_rs.csv", row.names = TRUE)
write.csv(R_ri, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R_ri.csv", row.names = TRUE)
write.csv(R_rs, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R_rs.csv", row.names = TRUE)



ri_results <- list(
  model_type = "random_intercept",
  fit_summary = subset(model_fit_summary, Model %in% c("m_ri")),
  parameter_summary = subset(parameter_summary, Model %in% c("m_ri")),
  treatment_all = list(
    D = D_ri,
    R = R_ri,
    fixed_summary = fixed_summary_ri,
    vcv_summary = vcv_summary_ri
  )
)

rs_results <- list(
  model_type = "random_intercept_random_slope",
  fit_summary = subset(model_fit_summary, Model %in% c("m_rs")),
  parameter_summary = subset(parameter_summary, Model %in% c("m_rs")),
  treatment_all = list(
    fit = m_rs,
    D = D_rs,
    R = R_rs,
    #fixed_summary = fixed_summary_rs,
    vcv_summary = vcv_summary_rs
  )
)

saveRDS(ri_results, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/random_intercept_results_num_cgi.rds")
saveRDS(rs_results, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/random_slope_results_num_cgi.rds")
end_time <- Sys.time()

nrow(m_ri$Sol)
nrow(m_rs$Sol)

plot(m_ri$Sol)
plot(m_rs$Sol)

library(coda)
effectiveSize(m_ri$Sol)
autocorr.diag(m_rs$Sol)

effectiveSize(m_rs$Sol)

effectiveSize(m_rs$VCV)
effectiveSize(m_rs$Sol)


##matrix_diagnostics_D
matrix_diagnostics_D <- function(M, model_name, matrix_name) {
  M <- as.matrix(M)
  M <- 0.5 * (M + t(M))  # numerical symmetry
  
  eigvals <- tryCatch(
    eigen(M, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(M))
  )
  
  det_val <- tryCatch(det(M), error = function(e) NA_real_)
  cond_val <- tryCatch(kappa(M, exact = TRUE), error = function(e) NA_real_)
  rcond_val <- if (is.na(cond_val)) NA_real_ else 1 / cond_val
  
  chol_ok <- !inherits(
    tryCatch(chol(M), error = function(e) e),
    "error"
  )
  
  status <- if (!chol_ok) {
    "not PD"
  } else if (!is.na(rcond_val) && rcond_val < 1e-8) {
    "PD but numerically singular"
  } else if (!is.na(rcond_val) && rcond_val < 1e-4) {
    "PD but ill-conditioned"
  } else {
    "PD"
  }
  
  data.frame(
    Model = model_name,
    Matrix = matrix_name,
    Determinant = det_val,
    Condition_number = cond_val,
    Status = status,
    stringsAsFactors = FALSE
  )
}

diag_table <- rbind(
  matrix_diagnostics_D(D_ri_out$blocks[[1]], "Endpoint-specific RI", "D0"),
  matrix_diagnostics_D(D_ri_out$blocks[[2]], "Endpoint-specific RI", "D1"),
  matrix_diagnostics_D(D_rs_out$blocks[[1]], "RI + RS", "D0"),
  matrix_diagnostics_D(D_rs_out$blocks[[2]], "RI + RS", "D1")
)

diag_table

