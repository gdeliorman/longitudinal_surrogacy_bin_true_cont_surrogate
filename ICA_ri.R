##ICA calculation with random intercept (no random slope)

library(Surrogate)
library(matrixcalc) 
library(fastmatrix)

##read outputs
ri_result<- readRDS(file= "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/random_intercept_results_num_cgi.rds" )
V_ri_matrix<- read.csv(file="/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R_ri.csv", row.names = 1, check.names = FALSE)
D_ri_matrix<- read.csv(file="/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/D_ri.csv", row.names = 1, check.names = FALSE)

D_ri_matrix[D_ri_matrix == 0] <- NA


V_ri_matrix <- as.matrix(V_ri_matrix)
D_ri_matrix <- as.matrix(D_ri_matrix)


storage.mode(V_ri_matrix) <- "numeric"
storage.mode(D_ri_matrix) <- "numeric"


D_ri_matrix_corr <- cov2cor(as.matrix(D_ri_matrix))

Sigma_1 <- diag(4)
Sigma_2<- V_ri_matrix
Sigma_3<- V_ri_matrix
Sigma_3[3, 4] <- Sigma_3[4, 3] <- NA



ICA_D<- ICA.ContCont(T0S0=D_ri_matrix_corr[1,3], T1S1=D_ri_matrix_corr[2,4], 
                     T0T0=D_ri_matrix[1,1], T1T1=D_ri_matrix[2,2], 
                     S0S0=D_ri_matrix[3,3], S1S1=D_ri_matrix[4,4], 
                     T0T1=seq(-1, 1, by=0.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))


ICA_D_pos_def<- ICA_D$Pos.Def
ICA_Sigma_3<- ICA.ContCont(T0S0=0, T1S1=0, 
                           T0T0=1, T1T1=1, 
                           S0S0=Sigma_3[3,3], S1S1=Sigma_3[4,4], 
                           T0T1=0, T0S1=0, T1S0=0, S0S1=seq(-1, 1, by=.005))

ICA_Sigma_3_pos_def<- ICA_Sigma_3$Pos.Def

# Q<- matrix(c(-1, 0,1,0,
#              0, -1,0,1), nrow=2, ncol=4)
p<- 5 ##time points
time_points = c(1, 2, 4, 6, 8)

build_Z_ri <- function(p) {
  ## Row order: T0(1:p), T1(1:p), S0(1:p), S1(1:p)
  Z <- matrix(0, nrow = 4 * p, ncol = 4)
  Z[1:p, 1] <- 1
  Z[(p + 1):(2 * p), 2] <- 1
  Z[(2 * p + 1):(3 * p), 3] <- 1
  Z[(3 * p + 1):(4 * p), 4] <- 1
  colnames(Z) <- c("T0", "T1", "S0", "S1")
  Z
}

Z<- build_Z_ri(5)

idx_T0 <- 1:p
idx_T1 <- (p + 1):(2 * p)
idx_S0 <- (2 * p + 1):(3 * p)
idx_S1 <- (3 * p + 1):(4 * p)


R2_1 <- c()
R2_2 <- c()
R2_3 <- c()

# i<-1
# j<-1

#for (i in 7436:nrow(ICA_D_pos_def)) {
for (i in 1:nrow(ICA_D_pos_def)) {
    
    
    D_ri_matrix[1,2]<- D_ri_matrix[2,1]<- cov_d_T0T1<-  ICA_D_pos_def[i,1]* sqrt(D_ri_matrix[1,1]*D_ri_matrix[2,2])
    D_ri_matrix[1,4]<- D_ri_matrix[4,1]<- cov_d_T0S1<-  ICA_D_pos_def[i,3]* sqrt(D_ri_matrix[1,1]*D_ri_matrix[4,4])
    D_ri_matrix[2,3]<- D_ri_matrix[3,2]<- cov_d_T1S0<-  ICA_D_pos_def[i,4]* sqrt(D_ri_matrix[2,2]*D_ri_matrix[3,3])
    D_ri_matrix[3,4]<- D_ri_matrix[4,3]<- cov_d_S0S1<-  ICA_D_pos_def[i,6]* sqrt(D_ri_matrix[3,3]*D_ri_matrix[4,4])
    
    
    R1 <- matrix(0, nrow = 4 * p, ncol = 4 * p)
    
    R1[c(1, 6, 11, 16), c(1, 6, 11, 16)] <- Sigma_1
    R1[c(2, 7, 12, 17), c(2, 7, 12, 17)] <- Sigma_1
    R1[c(3, 8, 13, 18), c(3, 8, 13, 18)] <- Sigma_1
    R1[c(4, 9, 14, 19), c(4, 9, 14, 19)] <- Sigma_1
    R1[c(5, 10, 15, 20), c(5, 10, 15, 20)] <- Sigma_1
    
    R2 <- matrix(0, nrow = 4 * p, ncol = 4 * p)
    
    R2[c(1, 6, 11, 16), c(1, 6, 11, 16)] <- Sigma_2
    R2[c(2, 7, 12, 17), c(2, 7, 12, 17)] <- Sigma_2
    R2[c(3, 8, 13, 18), c(3, 8, 13, 18)] <- Sigma_2
    R2[c(4, 9, 14, 19), c(4, 9, 14, 19)] <- Sigma_2
    R2[c(5, 10, 15, 20), c(5, 10, 15, 20)] <- Sigma_2
    
    Sigma_y_1 <- Z %*% D_ri_matrix %*% t(Z) + R1
    Sigma_y_2 <- Z %*% D_ri_matrix %*% t(Z) + R2
    
    
    Sigma_TT_1 <- Sigma_y_1[idx_T1, idx_T1] + Sigma_y_1[idx_T0, idx_T0] -
      Sigma_y_1[idx_T1, idx_T0] - Sigma_y_1[idx_T0, idx_T1]
    
    Sigma_SS_1 <- Sigma_y_1[idx_S1, idx_S1] + Sigma_y_1[idx_S0, idx_S0] -
      Sigma_y_1[idx_S1, idx_S0] - Sigma_y_1[idx_S0, idx_S1]
    
    Sigma_TS_1 <- Sigma_y_1[idx_T1, idx_S1] - Sigma_y_1[idx_T1, idx_S0] -
      Sigma_y_1[idx_T0, idx_S1] + Sigma_y_1[idx_T0, idx_S0]
    
    Sigma_Delta_1 <- rbind(
      cbind(Sigma_TT_1, Sigma_TS_1),
      cbind(t(Sigma_TS_1), Sigma_SS_1)
     )
    
    Lambda_1 <- det(Sigma_Delta_1) / (det(Sigma_TT_1) * det(Sigma_SS_1))
    R2_Lambda_1 <- 1 - Lambda_1
    #R2_1 <- c(R2_1, R2_Lambda_1)
    
    Sigma_TT_2 <- Sigma_y_2[idx_T1, idx_T1] + Sigma_y_2[idx_T0, idx_T0] -
      Sigma_y_2[idx_T1, idx_T0] - Sigma_y_2[idx_T0, idx_T1]
    
    Sigma_SS_2 <- Sigma_y_2[idx_S1, idx_S1] + Sigma_y_2[idx_S0, idx_S0] -
      Sigma_y_2[idx_S1, idx_S0] - Sigma_y_2[idx_S0, idx_S1]
    
    Sigma_TS_2 <- Sigma_y_2[idx_T1, idx_S1] - Sigma_y_2[idx_T1, idx_S0] -
      Sigma_y_2[idx_T0, idx_S1] + Sigma_y_2[idx_T0, idx_S0]
    
    Sigma_Delta_2 <- rbind(
      cbind(Sigma_TT_2, Sigma_TS_2),
      cbind(t(Sigma_TS_2), Sigma_SS_2)
    )
    
    Lambda_2 <- det(Sigma_Delta_2) / (det(Sigma_TT_2) * det(Sigma_SS_2))
    R2_Lambda_2<- 1 - Lambda_2
    #R2_2 <- c(R2_2, R2_Lambda_2)
    
    write(R2_Lambda_1, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_1.txt", append = TRUE)
    write(R2_Lambda_2, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_2.txt", append = TRUE)
    
    for (j in 1:nrow(ICA_Sigma_3_pos_def)) {
    
    Sigma_3[3,4] <- Sigma_3[4,3] <- cov_sigma_S0S1<-  ICA_Sigma_3_pos_def[j,6]* sqrt(Sigma_3[3,3]*Sigma_3[4,4])
  
    
    R3 <- matrix(0, nrow = 4 * p, ncol = 4 * p)
    
    R3[c(1, 6, 11, 16), c(1, 6, 11, 16)] <- Sigma_3
    R3[c(2, 7, 12, 17), c(2, 7, 12, 17)] <- Sigma_3
    R3[c(3, 8, 13, 18), c(3, 8, 13, 18)] <- Sigma_3
    R3[c(4, 9, 14, 19), c(4, 9, 14, 19)] <- Sigma_3
    R3[c(5, 10, 15, 20), c(5, 10, 15, 20)] <- Sigma_3
    
    Sigma_y_3 <- Z %*% D_ri_matrix %*% t(Z) + R3
    
    Sigma_TT_3 <- Sigma_y_3[idx_T1, idx_T1] + Sigma_y_3[idx_T0, idx_T0] -
      Sigma_y_3[idx_T1, idx_T0] - Sigma_y_3[idx_T0, idx_T1]
    
    Sigma_SS_3 <- Sigma_y_3[idx_S1, idx_S1] + Sigma_y_3[idx_S0, idx_S0] -
      Sigma_y_3[idx_S1, idx_S0] - Sigma_y_3[idx_S0, idx_S1]
    
    Sigma_TS_3 <- Sigma_y_3[idx_T1, idx_S1] - Sigma_y_3[idx_T1, idx_S0] -
      Sigma_y_3[idx_T0, idx_S1] + Sigma_y_3[idx_T0, idx_S0]
    
    Sigma_Delta_3 <- rbind(
      cbind(Sigma_TT_3, Sigma_TS_3),
      cbind(t(Sigma_TS_3), Sigma_SS_3)
    )
    
    
    Lambda_3 <- det(Sigma_Delta_3) / (det(Sigma_TT_3) * det(Sigma_SS_3))
    R2_Lambda_3 <- 1 - Lambda_3
    #R2_3 <- c(R2_3, R2_Lambda_3)
    
    write(R2_Lambda_3, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_3.txt", append = TRUE)
    
    
    
    }
}

R2_1 <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_1.txt",
             quiet = TRUE)

R2_2 <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_2.txt",
             quiet = TRUE)

R2_3 <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_3.txt",
             quiet = TRUE)

length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_1.txt"))
length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_2.txt"))
length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_ri_3.txt"))


summary_fun <- function(x) {
  c(
    n = length(x),
    min = min(x),
    Q1 = unname(quantile(x, 0.25)),
    Q2 = median(x),
    mean = mean(x),
    Q3 = unname(quantile(x, 0.75)),
    max = max(x),
    IQR = IQR(x)
  )
}

ri_summary <- rbind(
  Sigma_1 = summary_fun(R2_1),
  Sigma_2 = summary_fun(R2_2),
  Sigma_3 = summary_fun(R2_3)
)


boxplot(R2_1, R2_2, R2_3,
        names = expression(Sigma[1],  Sigma[2], Sigma[3]),
        ylim = c(0, 1),
        ylab = expression(R[Lambda]^2),
        col = "lightblue")


par(mfrow = c(1, 1))

hist(R2_1, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[1]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[1])),
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91", border = "white")

hist(R2_2, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[2]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[2])),
     
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91", border = "white")

hist(R2_3, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[3]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[3])),
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91",
     border = "white")
box()


