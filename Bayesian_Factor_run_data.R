

rm(list=ls(all=TRUE))
gc(reset = TRUE)

library(lpSolve)
library('data.table')
library(truncnorm)


source("Bayesian_DEA_MCMC_functions.R")


###############################################
###############################################


data <- data.table(read.table("mkt2015-data.txt", header=TRUE))
head(data)

myYear <- 2001
myYears <- 2001:2010

log_BF_all_years <- list()

for(myYear in myYears){
    cat(paste0("Current year is ",myYear,"\n"))
    XX <- as.matrix(data[data$year == myYear, paste0("x",1:5)])
    YY <- as.matrix(data[data$year == myYear, paste0("y",1:5)])
    WW <- as.matrix(data[data$year == myYear, paste0("w",1:5)])
    
    k = nrow(XX)
    
    P_I = WW
    
    m = ncol(YY)
    n = ncol(XX)
    
    setwd(paste0(getwd(),"/VRS/",myYear))

    TOTAL_log_BF = 0 
    BIG = 1.0e+15
    TOL = 1
    
    Method = 0
    BiM = 1.0e+6
    
    k1 = 210
    k1 = k
    
    log_BF_all <- matrix(NA, nrow = k1, ncol = 3)
    colnames(log_BF_all) <- c("Bank", "C_VRS", "NC_VRS")
    
    for(kkk in 1:k1)
    {
        gc(reset = TRUE)
        
        filename = paste('Results_posterior_', kkk, '.RData', sep = '')
        load(filename )
        
        Z_save = out_list[[1]]
        InnerRej = out_list[[2]] 
        MCMCRej = out_list[[3]] 
        
        if(length(Z_save) > k-1) 
        {	
            Z_save = Z_save[seq(1, nrow(Z_save), by = 10),]
            
            SIM = nrow(Z_save)
            
            
            #---------------------------------------------
            # Cost minimization:
            #---------------------------------------------
            
            c_k = c(P_I[kkk,]%*%t(XX))
            
            Max_c = max(abs(c_k))
            c_k = c_k/Max_c
            
            w_tilde_o <- c_k[kkk]
            w_tilde <- c_k[-kkk]
            
            c_k = c_k[-kkk] - w_tilde_o
            
            HIST = w_tilde_o + as.numeric((Z_save)%*%c_k)
            d = density(Max_c*HIST[HIST >= quantile(HIST, 0.15)], bw = Max_c/1000000)
            
            par(mfrow=c(1,2))
            hist(Max_c*HIST[HIST >= quantile(HIST, 0.15)], main = '', breaks = 50, xlab = 'Cost')
            plot(d, main = '', xlab = 'Cost')
            
            if(Method == 0) # Direct enumeration from the MCMC simulation
            {		
                CC = cov(Z_save)
                MM = min(1.0e-20, quantile(diag(CC)[which(diag(CC) > 0)], prob = 0.0))
                IND = which(diag(CC) >= MM)
                
                Z_save = Z_save[,IND]
                Z_mean = colMeans(Z_save)
                
                RowCheck = function(row){sum((row >= TOL*(floor(BIG*Z_mean))/BIG) & (row <= (2-TOL)*(ceiling(BIG*Z_mean))/BIG) ) >= k-1}
                selected_rows <- apply(Z_save, 1, RowCheck)
                PSUM = max(1, sum(selected_rows))
                
                PP_D = log(PSUM/nrow(Z_save))
                
                DENOM_log_exact = -as.numeric(c_k[IND]%*%Z_mean) - w_tilde_o + PP_D
                
            }
            if(Method == 1) # Gaussian approximation
            {			
                CC = cov(Z_save)
                MM = min(1.0e-20, quantile(diag(CC)[which(diag(CC) > 0)], prob = 0.10))
                IND = which(diag(CC) >= MM)
                CC_reduced = CC[IND, IND]        
                #V_det = det(CC_reduced)
                
                Z_save = Z_save[,IND]
                Z_mean = colMeans(Z_save)
                
                if(length(IND) == 0)
                {
                    PP_D = 0.5*(log(2*pi) + log(max(diag(CC))))
                }
                if(length(IND) == 1)
                {
                    PP_D = 0.5*(log(2*pi) + log(CC[1,1]))
                }else
                {
                    V_det_log <- determinant(CC_reduced)
                    PP_D = ((length(IND)/2)*log(2*pi)) + 0.5*unlist(V_det_log)[1]
                }
                
                if(PP_D <= -Inf)
                {
                    CC_reduced = CC_reduced + 1.1*diag(diag(CC_reduced))
                    
                    if(length(IND) == 0)
                    {
                        PP_D = 0.5*(log(2*pi) + log(max(diag(CC))))
                    }
                    if(length(IND) == 1)
                    {
                        PP_D = 0.5*(log(2*pi) + log(CC[1,1]))
                    }else
                    {
                        V_det_log <- determinant(CC_reduced)
                        PP_D = ((length(IND)/2)*log(2*pi)) + 0.5*unlist(V_det_log)[1]
                    }
                    
                }
                
                DENOM_log_exact = -as.numeric(c_k[IND]%*%Z_mean) - w_tilde_o + PP_D
                
            }
            if(Method == 2) # Dirichlet approximation
            {			
                Z_save_DIR = cbind(Z_save[,IND], 1 - rowSums(Z_save[,IND]))
                
                ALPHA = dirichlet.mle(Z_save_DIR)$alpha
                
                PP_D = sum(lfactorial(ALPHA-1)) - lfactorial(sum(ALPHA)-1)
                
                DENOM_log_exact = as.numeric(-c_k%*%Z_mean - w_tilde_o + PP_D)
                
            }
            if(Method == 3) # Multinomial approximation
            {			
                Z_save_MUL = round(Z_save * BiM)
                Z_save_MUL = cbind(Z_save_MUL, BiM - rowSums(Z_save_MUL))
                
                PPi = (colSums(Z_save_MUL)/nrow(Z_save_MUL))/BiM
                
                Z_mean_MUL = unlist(apply(Z_save_MUL, 2, median ))
                
                PP_D = dmultinom(x = Z_mean, size = sum(Z_mean_MUL), prob = PPi, log = TRUE)
                
                DENOM_log_exact = as.numeric(-c_k%*%Z_mean - w_tilde_o + PP_D)
                
            }
            
            #---------------------------------------------
            
            y_tilde <- numeric(k-1)
            for(h in 1:(k-1))
            {
                y_tilde[h] <- max( YY[kkk,] / YY[-kkk,][h,] )
            }
            NUM = exp(-w_tilde_o) + sum( exp(-w_tilde ) * (y_tilde < 1)  )
            
            currentK <- formatC(kkk, format = "d", width = 3)
            BF_exact_log <- log(NUM)- DENOM_log_exact
            BF_exact_log_ <- formatC(BF_exact_log, width = 6, format = "f", digits = 3)
            InnerRej_ <- formatC(InnerRej, format = "f", digits = 3)
            MCMCRej_ <- formatC(MCMCRej, format = "f", digits = 3)
            
            cat(paste0("Activity k = ",currentK, "; log_BF = ", BF_exact_log_, ", InnerRej = ", InnerRej_, ", MCMCRej = ", MCMCRej_,  "\n"))
            
            TOTAL_log_BF = TOTAL_log_BF + BF_exact_log
            
            # c("Bank", "C_VRS", "NC_VRS")
            log_BF_all[kkk,] <- c(kkk, DENOM_log_exact, log(NUM))
            print(log_BF_all[kkk,])
            
        } else {
            log_BF_all[kkk,] <- c(kkk, 0, 0)
        }
    }
    
    log_BF_all_years[[myYear]] <- log_BF_all
    
    TOTAL_log_BF
    
    print( mean( log_BF_all[,3] - log_BF_all[,2] ) )
    
}

rm(list=setdiff(ls(), "log_BF_all_years"))

save.image(file = "log_BF_all_years.RData")
