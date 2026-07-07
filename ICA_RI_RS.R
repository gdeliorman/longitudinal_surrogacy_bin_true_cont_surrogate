##ICA calculation with random effects (no random slope)

library(Surrogate)
library(matrixcalc) 
library(fastmatrix)

##read outputs

rs_result<- readRDS(file= "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/random_slope_results_num_cgi.rds" )
V_rs_matrix<- read.csv(file="/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R_rs.csv", row.names = 1, check.names = FALSE)
D_rs_matrix<- read.csv(file="/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/D_rs.csv", row.names = 1, check.names = FALSE)


V_rs_matrix <- as.matrix(V_rs_matrix)
D_rs_matrix <- as.matrix(D_rs_matrix)


storage.mode(V_rs_matrix) <- "numeric"
storage.mode(D_rs_matrix) <- "numeric"

D_rs_matrix[D_rs_matrix == 0] <- NA
D_rs_matrix_corr <- cov2cor(as.matrix(D_rs_matrix))

Sigma_1 <- diag(4)
Sigma_2<- V_rs_matrix
Sigma_3<- V_rs_matrix
Sigma_3[3, 4] <- Sigma_3[4, 3] <- NA

p<- 5 ##time points
time_points = c(1, 2, 4, 6, 8)

ICA_Sigma_3<- ICA.ContCont(T0S0=0, T1S1=0, 
                           T0T0=1, T1T1=1, 
                           S0S0=Sigma_3[3,3], S1S1=Sigma_3[4,4], 
                           T0T1=0, T0S1=0, T1S0=0, S0S1=seq(-1, 1, by=.005))

ICA_Sigma_3_pos_def<- ICA_Sigma_3$Pos.Def



#ica_mpc<-ICA.ContCont.MultS.MPC(M=10000,N=653,Sigma=D_rs_matrix,prob = NULL,Seed=123,
 #                               Save.Corr=TRUE, Show.Progress=FALSE)

ica_mpc<-ICA.ContCont.MultS.MPC(M=50000,N=568,Sigma=D_rs_matrix,prob = NULL,Seed=123,
                                Save.Corr=TRUE, Show.Progress=FALSE)

# ica_mpc<-ICA.ContCont.MultS.MPC(M=100000,N=568,Sigma=D_rs_matrix,prob = NULL,Seed=123,
#                                 Save.Corr=TRUE, Show.Progress=FALSE)
# 

valid_rows <- which(complete.cases(ica_mpc$Lower.Dig.Corrs.All))
complete_D<- list()

R_template<- D_rs_matrix_corr

for (i in valid_rows) {
  
  R<- R_template
  Sigma_D<- D_rs_matrix
  
  ##fill R matrix
  R[3,1]<- R[1,3]<-  ica_mpc$Lower.Dig.Corrs.All[i,2]
  R[4,1]<- R[1,4]<-  ica_mpc$Lower.Dig.Corrs.All[i,3]
  R[7,1]<- R[1,7]<-  ica_mpc$Lower.Dig.Corrs.All[i,6]
  R[8,1]<- R[1,8]<-  ica_mpc$Lower.Dig.Corrs.All[i,7]
  
  R[2,3]<- R[3,2]<-  ica_mpc$Lower.Dig.Corrs.All[i,8]
  R[2,4]<- R[4,2]<-  ica_mpc$Lower.Dig.Corrs.All[i,9]
  R[2,7]<- R[7,2]<-  ica_mpc$Lower.Dig.Corrs.All[i,12]
  R[2,8]<- R[8,2]<-  ica_mpc$Lower.Dig.Corrs.All[i,13]
  
  R[3,5]<- R[5,3]<-  ica_mpc$Lower.Dig.Corrs.All[i,15]
  R[3,6]<- R[6,3]<-  ica_mpc$Lower.Dig.Corrs.All[i,16]
  R[4,5]<- R[5,4]<-  ica_mpc$Lower.Dig.Corrs.All[i,19]
  R[4,6]<- R[6,4]<-  ica_mpc$Lower.Dig.Corrs.All[i,20]
  
  R[7,5]<- R[5,7]<-  ica_mpc$Lower.Dig.Corrs.All[i,24]
  R[8,5]<- R[5,8]<-  ica_mpc$Lower.Dig.Corrs.All[i,25]
  R[7,6]<- R[6,7]<-  ica_mpc$Lower.Dig.Corrs.All[i,26]
  R[8,6]<- R[6,8]<-  ica_mpc$Lower.Dig.Corrs.All[i,27]
  
  R <- (R + t(R)) / 2
  
  if( is.positive.definite(R) ==TRUE){
    Sigma_D[3,1]<- Sigma_D[1,3]<-  R[1,3]*sqrt(Sigma_D[1,1]*Sigma_D[3,3])
    Sigma_D[4,1]<- Sigma_D[1,4]<-  R[1,4]*sqrt(Sigma_D[1,1]*Sigma_D[4,4])
    Sigma_D[7,1]<- Sigma_D[1,7]<-  R[1,7]*sqrt(Sigma_D[1,1]*Sigma_D[7,7])
    Sigma_D[8,1]<- Sigma_D[1,8]<-  R[1,8]*sqrt(Sigma_D[1,1]*Sigma_D[8,8])
    
    Sigma_D[2,3]<- Sigma_D[3,2]<-  R[2,3]*sqrt(Sigma_D[2,2]*Sigma_D[3,3])
    Sigma_D[2,4]<- Sigma_D[4,2]<-  R[2,4]*sqrt(Sigma_D[2,2]*Sigma_D[4,4])
    Sigma_D[2,7]<- Sigma_D[7,2]<-  R[2,7]*sqrt(Sigma_D[2,2]*Sigma_D[7,7])
    Sigma_D[2,8]<- Sigma_D[8,2]<-  R[2,8]*sqrt(Sigma_D[2,2]*Sigma_D[8,8])
    
    Sigma_D[3,5]<- Sigma_D[5,3]<-  R[3,5]*sqrt(Sigma_D[5,5]*Sigma_D[3,3])
    Sigma_D[3,6]<- Sigma_D[6,3]<-  R[3,6]*sqrt(Sigma_D[6,6]*Sigma_D[3,3])
    Sigma_D[4,5]<- Sigma_D[5,4]<-  R[4,5]*sqrt(Sigma_D[5,5]*Sigma_D[4,4])
    Sigma_D[4,6]<- Sigma_D[6,4]<-  R[4,6]*sqrt(Sigma_D[6,6]*Sigma_D[4,4])
    
    Sigma_D[7,5]<- Sigma_D[5,7]<-  R[7,5]*sqrt(Sigma_D[5,5]*Sigma_D[7,7])
    Sigma_D[8,5]<- Sigma_D[5,8]<-  R[8,5]*sqrt(Sigma_D[5,5]*Sigma_D[8,8])
    Sigma_D[7,6]<- Sigma_D[6,7]<-  R[7,6]*sqrt(Sigma_D[6,6]*Sigma_D[7,7])
    Sigma_D[8,6]<- Sigma_D[6,8]<-  R[8,6]*sqrt(Sigma_D[6,6]*Sigma_D[8,8])
    complete_D[[i]]<- Sigma_D
  }
}

complete_D<-complete_D[!sapply(complete_D,is.null)]
##check D are unique?
length(unique(complete_D))
length(complete_D)



# Q<- matrix(c(-1, 0,1,0,
#              0, -1,0,1), nrow=2, ncol=4)


Z <- matrix(0, nrow = 4 * p, ncol = 8)

## T0 rows
Z[1:p, 1] <- 1
Z[1:p, 2] <- time_points

## T1 rows
Z[(p + 1):(2 * p), 3] <- 1
Z[(p + 1):(2 * p), 4] <- time_points

## S0 rows
Z[(2 * p + 1):(3 * p), 5] <- 1
Z[(2 * p + 1):(3 * p), 6] <- time_points

## S1 rows
Z[(3 * p + 1):(4 * p), 7] <- 1
Z[(3 * p + 1):(4 * p), 8] <- time_points

colnames(Z) <- c(
  "T0_int", "T0_slope",
  "T1_int", "T1_slope",
  "S0_int", "S0_slope",
  "S1_int", "S1_slope"
)


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

  
for (i in 1:length(complete_D)) {
  
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
  
  Sigma_y_1 <- Z %*% complete_D[[i]] %*% t(Z) + R1
  Sigma_y_2 <- Z %*% complete_D[[i]] %*% t(Z) + R2
  
  
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

  write(R2_Lambda_1, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_1.txt", append = TRUE)
  write(R2_Lambda_2, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_2.txt", append = TRUE)
  
  for (j in 1:nrow(ICA_Sigma_3_pos_def)) {
    
    Sigma_3[3,4] <- Sigma_3[4,3] <- cov_sigma_S0S1<-  ICA_Sigma_3_pos_def[j,6]* sqrt(Sigma_3[3,3]*Sigma_3[4,4])
    
    
    R3 <- matrix(0, nrow = 4 * p, ncol = 4 * p)
    
    R3[c(1, 6, 11, 16), c(1, 6, 11, 16)] <- Sigma_3
    R3[c(2, 7, 12, 17), c(2, 7, 12, 17)] <- Sigma_3
    R3[c(3, 8, 13, 18), c(3, 8, 13, 18)] <- Sigma_3
    R3[c(4, 9, 14, 19), c(4, 9, 14, 19)] <- Sigma_3
    R3[c(5, 10, 15, 20), c(5, 10, 15, 20)] <- Sigma_3
    
    Sigma_y_3 <- Z %*% complete_D[[i]] %*% t(Z) + R3
    
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

    write(R2_Lambda_3, "/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_3.txt", append = TRUE)
    
    
    
  }
}

R2_1_rs <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_1.txt",
             quiet = TRUE)

R2_2_rs <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_2.txt",
             quiet = TRUE)

R2_3_rs <- scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_3.txt",
             quiet = TRUE)

length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_1.txt"))
length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_2.txt"))
length(scan("/Users/gokcedeliorman/Desktop/new_study_phd/joint_model_output_2_new/R2_rs_3.txt"))


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
  Sigma_1 = summary_fun(R2_1_rs),
  Sigma_2 = summary_fun(R2_2_rs),
  Sigma_3 = summary_fun(R2_3_rs)
)

ri_summary


boxplot(R2_1_rs, R2_2_rs, R2_3_rs,
        names = expression(Sigma[1],  Sigma[2], Sigma[3]),
        ylim = c(0, 1),
        ylab = expression(R[Lambda]^2),
        col = "lightblue")


par(mfrow = c(1, 1))

hist(R2_1_rs, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[1]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[1])),
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91", border = "white")

hist(R2_2_rs, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[2]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[2])),
     
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91", border = "white")

hist(R2_3_rs, breaks = 40, xlim = c(0, 1),
     #main = expression(Sigma[3]),
     main = expression(paste(R[Lambda]^2, " for ", Sigma[3])),
     xlab = expression(R[Lambda]^2),
     col = "#7F8A91",
     border = "white")
box()


