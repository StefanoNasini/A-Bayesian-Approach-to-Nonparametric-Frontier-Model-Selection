

library(lpSolve)
library(truncnorm)


###############################################
###############################################


#---------------------------------------------
# Functions
#---------------------------------------------
# Generalized data structure for Cost min and Rev maxf
#---------------------------------------------

A_mat = function(s){

	if(s >= 1){
		out = t(YY)
	}else{
		out = -t(XX)
	}

	return(out)

}

q_mat = function(s){

	if(s >= 1){
		out = m
	}else{
		out = n
	}

	return(out)

}

p_mat = function(s,kkk){

	if(s >= 1){
		out = as.numeric(P_I[kkk,])
	}else{
		out = -as.numeric(P_O[kkk,])
	}

	return(out)

}

D_mat = function(s){

	if(s >= 1){
		out = diag(n)
	}else{
		out = -diag(m)
	}

	return(out)

}


b_mat = function(s, kkk){

	if(s >= 1){
		out = as.numeric(YY[kkk,])
	}else{
		out = -as.numeric(XX[kkk,])
	}

	return(out)

}



LP_mat = function(s, MM, dd){

	out = rbind( 
		cbind(A_mat(s)*(dd/MM), matrix(0, q_mat(s), q_mat(1-s))),
		cbind(A_mat(1-s)*(dd/MM), D_mat(s)),
		c(rep(1, k), rep(0, q_mat(1-s))),
		cbind(diag(rep(1, k)), matrix(0, k, q_mat(1-s)))
	)

	return(out)

}

LP_mat_eq = function(s, MM, dd, slack){

	if(slack){

		out = rbind( 
			cbind(A_mat(s)*(dd/MM), -diag(q_mat(s))),
			c(rep(1, k), rep(0, q_mat(s)))
		)

	}else{

		out = rbind( 
			cbind(A_mat(s)*(dd/MM)),
			c(rep(1, k))
		)

	}

	return(out)

}


LP_mat_eq_delta = function(s, MM, slack){

	if(slack){

		out = rbind( 
			cbind(A_mat(s)*(1/MM), -diag(q_mat(s)), rep(0, q_mat(s) ) ),
			c(rep(1, k), rep(0, q_mat(s)), -1 )
		)

	}else{

		out = rbind( 
			cbind(A_mat(s)*(1/MM), rep(0, q_mat(s) ) ),
			c(rep(1, k), -1)
		)

	}

	return(out)
	
}



g = function(zz_next, zz_old, ss){

	prob_candidate = 0
	dimension = length(zz_old)
	
	uub = cumsum(zz_next)[1:dimension]

	for(h in 1:dimension)
	{
		if(h <= 1)
		{
			prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a=0, b = 1, mean = zz_old[h], sd = ss))
		
		}else
		{

			if(uub[h] < 1 - 1.0e-10)
			{
				prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a=0, b = 1 - uub[h-1], mean = zz_old[h], sd = ss))
		
			}	
		}	


	}

	return(prob_candidate) 
}



g_delta_CRS = function(zz_next, zz_old, delta_next, delta_old, ss)
{
	dimension = length(zz_old)
	
	uub = cumsum(zz_next)[1:dimension]

	prob_candidate = log(dtruncnorm(delta_next, a=0, b = Inf, mean = delta_old, sd = sqrt(dimension)*ss))

	for(h in 1:dimension)
	{
		if(h <= 1)
		{
			prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next, mean = zz_old[h], sd = ss))
		}
		else
		{
			if(uub[h] < delta_next - 1.0e-10)
			{
				prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next - uub[h-1], mean = zz_old[h], sd = ss))

			}	
		}	
	}

	return(prob_candidate) 
}


g_delta_NIRS = function(zz_next, zz_old, delta_next, delta_old, ss)
{
	dimension = length(zz_old)
	
	uub = cumsum(zz_next)[1:dimension]

	prob_candidate = log(dtruncnorm(delta_next, a=0, b = 1, mean = delta_old, sd = sqrt(dimension)*ss))

	for(h in 1:dimension)
	{
		if(h <= 1)
		{
			prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next, mean = zz_old[h], sd = ss))
		}
		else
		{
			if(uub[h] < delta_next - 1.0e-10)
			{
				prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next - uub[h-1], mean = zz_old[h], sd = ss))
		
			}	
		}	
	}

	return(prob_candidate) 
}


g_delta_NDRS = function(zz_next, zz_old, delta_next, delta_old, ss)
{
	dimension = length(zz_old)
	
	uub = cumsum(zz_next)[1:dimension]

	prob_candidate = log(dtruncnorm(delta_next, a=1, b = Inf, mean = delta_old, sd = sqrt(dimension)*ss))

	for(h in 1:dimension)
	{
		if(h <= 1)
		{
			prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next, mean = zz_old[h], sd = ss))
		}
		else
		{
			if(uub[h] < delta_next - 1.0e-15)
			{
				prob_candidate = prob_candidate + log(dtruncnorm(zz_next[h], a = 0, b = delta_next - uub[h-1], mean = zz_old[h], sd = ss))
		
			}	
		}	
	}

	return(prob_candidate) 
}


Simulate_z = function(zz_old, zz_feasible, ss, A_kk, H_set, b_k){

	Reject = 0
	Count_reject = 0 
	out = list()

	dimension = length(zz_old)
	zz_next = rep(0,dimension)

	while(Reject == 0)
	{

		ub = 0

		for(h in 1:dimension)
		{
			if(ub < 1){
				zz_next[h] = rtruncnorm(1, a=0, b=(1 - ub), mean = zz_old[h], sd = ss)
				ub = ub + zz_next[h]
			}else{
				zz_next[h] = 0
			}

		}

		Reject = min(A_kk[, H_set]%*%zz_next >= b_k)

		Count_reject = Count_reject + 1
		
		if(Count_reject >= 100)
		{
			aa = runif(1, min=0, max=1)
			
			for(h in 1:dimension)
			{
				zz_next[h] = aa*zz_old[h] + (1-aa)*zz_feasible[h]
			}
			
			break
		}
	}

	out[[1]] = zz_next
	out[[2]] = Count_reject

	return(out) 
}


Simulate_z_CRS = function(zz_old, zz_feasible, delta_old, delta_feasible, ss, A_kk, H_set, b_k)
{
	Reject = 0
	Count_reject = 0 
	out = list()

	dimension = length(zz_old)
	zz_next = rep(0,dimension)
	
	while(Reject == 0)
	{
		delta_next = rtruncnorm(1, a=0, b=Inf, mean = delta_old, sd = sqrt(dimension)*ss)			
		ub = 0

		for(h in 1:dimension)
		{

			if(ub < delta_next)
			{
				zz_next[h] = rtruncnorm(1, a=0, b = max(1.0e-14/h, delta_next - ub), mean = zz_old[h], sd = ss)
				ub = ub + zz_next[h]
			}else{
				zz_next[h] = 0
			}

		}
		
		delta_next = max(delta_next, sum(zz_next))
		

		Reject = min(A_kk[, H_set]%*%zz_next >= (1-delta_next)*b_k)

		Count_reject = Count_reject + 1
		
		if(Count_reject >= 1000)
		{
			aa = runif(1, min=0, max=1)
			
			for(h in 1:dimension)
			{
				zz_next[h] = aa*zz_old[h] + (1-aa)*zz_feasible[h]
			}
			delta_next = aa*delta_old + (1-aa)*delta_feasible
			
			break
		}
	}

	out[[1]] = zz_next
	out[[2]] = delta_next
	out[[3]] = Count_reject

	return(out) 
}


Simulate_z_NIRS = function(zz_old, zz_feasible, delta_old, delta_feasible, ss, A_kk, H_set, b_k)
{
	Reject = 0
	Count_reject = 0 
	out = list()

	dimension = length(zz_old)
	zz_next = rep(0,dimension)
	
	while(Reject == 0)
	{
		delta_next = rtruncnorm(1, a=0, b=1, mean = delta_old, sd = sqrt(dimension)*ss)			
		ub = 0

		for(h in 1:dimension)
		{
			if(ub < delta_next)
			{
				zz_next[h] = rtruncnorm(1, a=0, b = max(1.0e-14/h, delta_next - ub), mean = zz_old[h], sd = ss)
				ub = ub + zz_next[h]
			}else{
				zz_next[h] = 0
			}

		}
		
		delta_next = max(delta_next, sum(zz_next))
		

		Reject = min(A_kk[, H_set]%*%zz_next >= (1-delta_next)*b_k)

		#print(cbind(A_kk[, H_set]%*%zz_next, b_k) )
		Count_reject = Count_reject + 1
		
		if(Count_reject >= 1000)
		{
			aa = runif(1, min=0, max=1)
			
			for(h in 1:dimension)
			{
				zz_next[h] = aa*zz_old[h] + (1-aa)*zz_feasible[h]
			}
			delta_next = aa*delta_old + (1-aa)*delta_feasible
			
			break
		}
	}

	out[[1]] = zz_next
	out[[2]] = delta_next
	out[[3]] = Count_reject

	return(out) 
}


Simulate_z_NDRS = function(zz_old, zz_feasible, delta_old, delta_feasible, ss, A_kk, H_set, b_k)
{
	Reject = 0
	Count_reject = 0 
	out = list()

	dimension = length(zz_old)
	zz_next = rep(0,dimension)
	
	while(Reject == 0)
	{
		delta_next = rtruncnorm(1, a=1, b=Inf, mean = delta_old, sd = sqrt(dimension)*ss)			
		ub = 0

		for(h in 1:dimension)
		{
			if(ub < delta_next)
			{
				zz_next[h] = rtruncnorm(1, a=0, b = max(1.0e-14/h, delta_next - ub), mean = zz_old[h], sd = ss)
				ub = ub + zz_next[h]
			}else{
				zz_next[h] = 0
			}

		}
		
		delta_next = max(delta_next, sum(zz_next))
		
		Reject = min(A_kk[, H_set]%*%zz_next >= (1-delta_next)*b_k)

		Count_reject = Count_reject + 1
		
		if(Count_reject >= 1000)
		{
			aa = runif(1, min=0, max=1)
			
			for(h in 1:dimension)
			{
				zz_next[h] = aa*zz_old[h] + (1-aa)*zz_feasible[h]
			}
			delta_next = aa*delta_old + (1-aa)*delta_feasible
			
			break
		}
	}

	out[[1]] = zz_next
	out[[2]] = max(1, delta_next)
	out[[3]] = Count_reject

	return(out) 
}


###########################################
# MH algorithm
###########################################

MetropolisHansongs_z_posterior = function(kk = 2, tolerance = 1.0e-15, ITER_max = 50000, ITER_min = 200, proposal_var = 0.01, technology = 0, H = 1000000000)
{

	MH_count_reject = 0
	Inner_count_reject = 0
	Z_list = matrix(0,ITER_max, k-1) 

	VAR_MH_old = 1000

	iter = 1
	GAP = 1000
		
	out_list = list()
	
	#---------------------------------------------
	# Feasible point:
	#---------------------------------------------

	A_kk = LP_mat_eq(1, 1, 1, slack = FALSE) 

	b_k = c(rep(0,m), 1)

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c_k[-kk] - c_kk

	c_k = max(c_k) - 0.9*c_k
	
	A_kk = A_kk[1:m,-kk] - A_kk[1:m,kk]

	AA_kk = rbind(A_kk, rep(1,k-1))

	DIR = c(rep(">=", m), "<=")

	out_EDM = lp (direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_feasible = out_EDM$solution[1:(k-1)]

	#---------------------------------------------
	# Cost minimization:
	#---------------------------------------------

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c_k[-kk] - c_kk

	out_EDM = lp (direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_opt = out_EDM$solution[1:(k-1)]
	
	#################################################################
	#
	# Obtaining the expected posterior by MCMC
	#
	#################################################################

	set.seed(1234)

	# Sort activity vector in decreasing order
	
	if(max(z_opt) > 0){		
		H_set= order(z_opt, decreasing = TRUE)
	}else{		
		H_set= 1:(k-1)			
	}
	
	if(sum(abs(z_feasible - z_opt)) <= 0)
	{
		Z_list = matrix(rep(z_feasible, each = ITER_max), nrow = ITER_max, byrow = TRUE)
		
		z_plot = z_feasible[1]

		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)
		
		plot(0, pch = '')
		legend("topleft", legend = c(paste("Expectation = ", z_plot, sep = ''), paste("Optimal solution = ", z_opt[1], sep = '')), col=c(3, 2), lwd=3, cex=1)

		dev.off()
			
		Z_save = z_opt
		
		out_list[[1]] = Z_save
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)	
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
	else
	{		
		z_feasible = (0.99*z_opt[H_set] + 0.01*z_feasible[H_set])
		z_opt = z_opt[H_set]
			
		# Initialization

		z_old = z_feasible

		sss = proposal_var

		while(((GAP >= tolerance) & (iter <= ITER_max)) || (iter < ITER_min) )
		{
			if((iter > 100) & ((iter/10) == round(iter/10)) )
			{			
				Rej_rate = MH_count_reject/(1 + iter)
				
				z_old = colMeans(Z_list[1:iter, ])
								
				if(Rej_rate > 0.5)
				{ 
					sss = max(1.0e-10, sss*(1 - Rej_rate*0.01))
				}
				else
				{
					sss = max(1.0e-10, sss*(1 + (1-Rej_rate)*0.01))
				}
			}

			Simul = Simulate_z(z_old, z_opt, sss, A_kk, H_set, rep(0,m))
			z_next = Simul[[1]]
			Count_reject = Simul[[2]]
			Inner_count_reject = Inner_count_reject +  Count_reject

			if(Count_reject >= 100)
			{
				prob_candidate_next = 1
				prob_candidate_old = 1
			}
			else
			{
				prob_candidate_next = g(z_next, z_old, sss) 
				prob_candidate_old = g(z_old, z_next, sss) 
			}
			
			AA = -H*as.numeric(c_k[H_set]%*%z_next) + H*as.numeric(c_k[H_set]%*%z_old) + prob_candidate_old - prob_candidate_next

			uu = runif(1)

			if(log(uu) <= AA)
			{
				z_old = z_next
			}
			else
			{
				MH_count_reject = MH_count_reject + 1
			}

			Z_list[iter, ] = z_old
			
			if((iter > MH_count_reject + 10) & ((iter/500) == round(iter/500)))
			{
				VAR_MH = sqrt(diag(cov(Z_list[1:iter, ])))

				GAP = mean(abs(VAR_MH - VAR_MH_old))/(1 + max(VAR_MH_old))

				VAR_MH_old = VAR_MH

			}

			cat(c(round(iter), round(Inner_count_reject/(1 + iter), 3),  round(MH_count_reject/(1 + iter), 3), round(mean(Z_list[1:iter, 1]),3) , max(VAR_MH_old), GAP, '\n'), sep = "\t") 

			iter = iter + 1

		}

		ITER_max = iter -1
		
		BURN_IN = round(ITER_max/20)
		BURN_IN = 1
	
		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)

		MCMC_ind = seq(BURN_IN, ITER_max, by = 1)
		z_plot = Z_list[MCMC_ind, 1]
		Z_save = Z_list[MCMC_ind, ]
		
		z_down = min(z_plot)
		z_up = max(z_plot)
			
		hist(z_plot, breaks = 100, ylab = "Probability", xlab = "z[1]", xlim = c(z_down, z_up), main = "", axes = FALSE)
		legend("topleft", legend=c(paste("Expectation = ", mean(z_plot), sep = ''), paste("Optimal solution = ", z_opt[1], sep = ''), paste("Rejection rate = ", round(MH_count_reject/(1 + iter), 3), sep = '') ), col=c(3, 2, 0), lwd=3, cex=1)
		
		dev.off()		
		
		out_list[[1]] = Z_save
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
			
	return(out_list)


}



MetropolisHansongs_z_posterior_CRS = function(kk = 2, tolerance = 1.0e-15, ITER_max = 50000, ITER_min = 200, proposal_var = 0.01, technology = 0, H = 1000000000)
{

	MH_count_reject = 0
	Inner_count_reject = 0
	Z_list = matrix(0,ITER_max, k) 

	VAR_MH_old = 1000

	iter = 1
	GAP = 1000
		
	out_list = list()
	
	#---------------------------------------------
	# Feasible point:
	#---------------------------------------------

	A_kk = A_mat(1) 

	b_k = c(A_kk[1:m,kk], 0)

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c_k[-kk] - c_kk

	c_k = c(max(c_k) - 0.9*c_k, 0)

	A_kk = cbind(A_kk[1:m,-kk] - A_kk[1:m,kk],  -A_kk[1:m,kk])
	
	AA_kk = rbind(A_kk, c(rep(1,k-1), -1))

	DIR = c(rep(">=", m), "<=")

	out_EDM = lp(direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_feasible = out_EDM$solution[1:(k-1)]
	delta_feasible = out_EDM$solution[k]
	
	#---------------------------------------------
	# Cost minimization:
	#---------------------------------------------

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c(c_k[-kk] - c_kk, c_kk)

	out_EDM = lp (direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_opt = out_EDM$solution[1:(k-1)]
	delta_opt = out_EDM$solution[k]
		
	#################################################################
	#
	# Obtaining the expected posterior by MCMC
	#
	#################################################################

	set.seed(1234)

	# Sort activity vector in decreasing order
	
	if(max(z_opt) > 0){		
		H_set= order(z_opt, decreasing = TRUE)
	}else{		
		H_set= 1:(k-1)			
	}
	
	if(sum(abs(z_feasible - z_opt)) <= 0)
	{
		Z_list = matrix(rep(z_feasible, each = ITER_max), nrow = ITER_max, byrow = TRUE)
		
		z_plot = z_feasible[1]

		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)
		
		plot(0, pch = '')
		legend("topleft", legend = c(paste("Expectation = ", z_plot, sep = ''), paste("Optimal solution = ", z_opt[1], sep = '')), col=c(3, 2), lwd=3, cex=1)

		dev.off()
		
		out_list[[1]] = c(z_opt, delta_opt[1])		
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)	
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
	else
	{		
		z_feasible = (0.99*z_opt[H_set] + 0.01*z_feasible[H_set])
		delta_feasible = (0.99*delta_opt + 0.01*delta_feasible)
		z_opt = z_opt[H_set]
			
		# Initialization

		z_old = z_feasible
		delta_old = delta_feasible
		
		sss = proposal_var

		while(((GAP >= tolerance) & (iter <= ITER_max)) || (iter < ITER_min) )
		{
			if((iter > 100) & ((iter/10) == round(iter/10)) )
			{			
				Rej_rate = MH_count_reject/(1 + iter)
				
				z_old = colMeans(Z_list[1:iter, (1:k-1)])
				delta_old = mean(Z_list[1:iter, k])
												
				if(Rej_rate > 0.5)
				{ 
					sss = max(1.0e-10, sss*(1 - Rej_rate*0.01))
				}
				else
				{
					sss = max(1.0e-10, sss*(1 + (1-Rej_rate)*0.01))
				}
			}

			Simul = Simulate_z_CRS(z_old, z_opt, delta_old, delta_opt, sss, A_kk, H_set, A_kk[1:m,kk])
			z_next = Simul[[1]]
			delta_next = Simul[[2]]
			Count_reject = Simul[[3]]
			
			Inner_count_reject = Inner_count_reject +  Count_reject

			if(Count_reject >= 1000)
			{
				prob_candidate_next = 1
				prob_candidate_old = 1
			}
			else
			{
				prob_candidate_next = g_delta_CRS(z_next, z_old, delta_next, delta_old, sss) 
				prob_candidate_old = g_delta_CRS(z_old, z_next, delta_old, delta_next, sss) 
			}
			
			AA = -H*as.numeric(c_k[H_set]%*%z_next + c_kk*delta_next) + H*as.numeric(c_k[H_set]%*%z_old + c_kk*delta_old) + prob_candidate_old - prob_candidate_next
			
			uu = runif(1)

			if(log(uu) <= AA)
			{
				z_old = z_next
				delta_old = delta_next
			}
			else
			{
				MH_count_reject = MH_count_reject + 1
			}

			Z_list[iter, ] = c(z_old, delta_old)
			
			if((iter > MH_count_reject + 10) & ((iter/500) == round(iter/500)))
			{
				VAR_MH = sqrt(diag(cov(Z_list[1:iter, ])))

				GAP = mean(abs(VAR_MH - VAR_MH_old))/(1 + max(VAR_MH_old))

				VAR_MH_old = VAR_MH

			}

			cat(c(round(iter), round(Inner_count_reject/(1 + iter), 3),  round(MH_count_reject/(1 + iter), 3), round(mean(Z_list[1:iter, 1]),3) , max(VAR_MH_old), GAP, '\n'), sep = "\t") 

			iter = iter + 1

		}

		ITER_max = iter -1
		
		BURN_IN = round(ITER_max/20)
		BURN_IN = 1
	
		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)

		MCMC_ind = seq(BURN_IN, ITER_max, by = 1)
		z_plot = Z_list[MCMC_ind, 1]
		Z_save = Z_list[MCMC_ind, ]
		
		z_down = min(z_plot)
		z_up = max(z_plot)
			
		hist(z_plot, breaks = 100, ylab = "Probability", xlab = "z[1]", xlim = c(z_down, z_up), main = "", axes = FALSE)
		legend("topleft", legend=c(paste("Expectation = ", mean(z_plot), sep = ''), paste("Optimal solution = ", z_opt[1], sep = ''), paste("Rejection rate = ", round(MH_count_reject/(1 + iter), 3), sep = '') ), col=c(3, 2, 0), lwd=3, cex=1)
		
		dev.off()		
		
		out_list[[1]] = Z_save
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
			
	return(out_list)


}



MetropolisHansongs_z_posterior_NIRS = function(kk = 2, tolerance = 1.0e-15, ITER_max = 50000, ITER_min = 200, proposal_var = 0.01, technology = 0, H = 1000000000)
{

	MH_count_reject = 0
	Inner_count_reject = 0
	Z_list = matrix(0,ITER_max, k) 

	VAR_MH_old = 1000

	iter = 1
	GAP = 1000
		
	out_list = list()
	
	#---------------------------------------------
	# Feasible point:
	#---------------------------------------------

	A_kk = A_mat(1) 

	b_k = c(A_kk[1:m,kk], 0, 1)

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c_k[-kk] - c_kk

	c_k = c(max(c_k) - 0.9*c_k, 0)

	A_kk = cbind(A_kk[1:m,-kk] - A_kk[1:m,kk],  -A_kk[1:m,kk])
	
	AA_kk = rbind(A_kk, c(rep(1,k-1), -1), c(rep(0,k-1), 1))

	DIR = c(rep(">=", m), "<=", "<=")

	out_EDM = lp(direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_feasible = out_EDM$solution[1:(k-1)]
	delta_feasible = out_EDM$solution[k]
	
	#---------------------------------------------
	# Cost minimization:
	#---------------------------------------------

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c(c_k[-kk] - c_kk, c_kk)

	out_EDM = lp (direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_opt = out_EDM$solution[1:(k-1)]
	delta_opt = out_EDM$solution[k]
		
	#################################################################
	#
	# Obtaining the expected posterior by MCMC
	#
	#################################################################

	set.seed(1234)

	# Sort activity vector in decreasing order
	
	if(max(z_opt) > 0){		
		H_set= order(z_opt, decreasing = TRUE)
	}else{		
		H_set= 1:(k-1)			
	}
	
	if(sum(abs(z_feasible - z_opt)) <= 0)
	{
		Z_list = matrix(rep(z_feasible, each = ITER_max), nrow = ITER_max, byrow = TRUE)
		
		z_plot = z_feasible[1]

		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)
		
		plot(0, pch = '')
		legend("topleft", legend = c(paste("Expectation = ", z_plot, sep = ''), paste("Optimal solution = ", z_opt[1], sep = '')), col=c(3, 2), lwd=3, cex=1)

		dev.off()
		
		out_list[[1]] = c(z_opt, delta_opt[1])		
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)	
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
	else
	{		
		delta_feasible = (0.99*delta_opt + 0.01*delta_feasible)
		z_opt = z_opt[H_set]
			
		# Initialization

		z_old = z_feasible
		delta_old = delta_feasible
		
		sss = proposal_var

		while(((GAP >= tolerance) & (iter <= ITER_max)) || (iter < ITER_min) )
		{
			if((iter > 100) & ((iter/10) == round(iter/10)) )
			{			
				Rej_rate = MH_count_reject/(1 + iter)
				
				z_old = colMeans(Z_list[1:iter, (1:k-1)])
				delta_old = mean(Z_list[1:iter, k])
												
				if(Rej_rate > 0.5)
				{ 
					sss = max(1.0e-10, sss*(1 - Rej_rate*0.01))
				}
				else
				{
					sss = max(1.0e-10, sss*(1 + (1-Rej_rate)*0.01))
				}
			}

			Simul = Simulate_z_NIRS(z_old, z_opt, delta_old, delta_opt, sss, A_kk, H_set, A_kk[1:m,kk])
			z_next = Simul[[1]]
			delta_next = Simul[[2]]
			Count_reject = Simul[[3]]
			
			Inner_count_reject = Inner_count_reject +  Count_reject

			if(Count_reject >= 1000)
			{
				prob_candidate_next = 1
				prob_candidate_old = 1
			}
			else
			{
				prob_candidate_next = g_delta_NIRS(z_next, z_old, delta_next, delta_old, sss) 
				prob_candidate_old = g_delta_NIRS(z_old, z_next, delta_old, delta_next, sss) 
			}
			
			AA = -H*as.numeric(c_k[H_set]%*%z_next + c_kk*delta_next) + H*as.numeric(c_k[H_set]%*%z_old + c_kk*delta_old) + prob_candidate_old - prob_candidate_next
			
			uu = runif(1)

			if(log(uu) <= AA)
			{
				z_old = z_next
				delta_old = delta_next
			}
			else
			{
				MH_count_reject = MH_count_reject + 1
			}

			Z_list[iter, ] = c(z_old, delta_old)
			
			#print(iter)
			#print(Z_list[iter, ])

			if((iter > MH_count_reject + 10) & ((iter/500) == round(iter/500)))
			{
				VAR_MH = sqrt(diag(cov(Z_list[1:iter, ])))

				GAP = mean(abs(VAR_MH - VAR_MH_old))/(1 + max(VAR_MH_old))

				VAR_MH_old = VAR_MH

			}

			cat(c(round(iter), round(Inner_count_reject/(1 + iter), 3),  round(MH_count_reject/(1 + iter), 3), round(mean(Z_list[1:iter, 1]),3) , max(VAR_MH_old), GAP, '\n'), sep = "\t") 

			iter = iter + 1

		}

		ITER_max = iter -1
		
		BURN_IN = round(ITER_max/20)
		BURN_IN = 1
	
		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)

		MCMC_ind = seq(BURN_IN, ITER_max, by = 1)
		z_plot = Z_list[MCMC_ind, 1]
		Z_save = Z_list[MCMC_ind, ]
		
		z_down = min(z_plot)
		z_up = max(z_plot)
			
		hist(z_plot, breaks = 100, ylab = "Probability", xlab = "z[1]", xlim = c(z_down, z_up), main = "", axes = FALSE)
		legend("topleft", legend=c(paste("Expectation = ", mean(z_plot), sep = ''), paste("Optimal solution = ", z_opt[1], sep = ''), paste("Rejection rate = ", round(MH_count_reject/(1 + iter), 3), sep = '') ), col=c(3, 2, 0), lwd=3, cex=1)
		
		dev.off()		
		
		out_list[[1]] = Z_save
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
			
	return(out_list)


}



MetropolisHansongs_z_posterior_NDRS = function(kk = 2, tolerance = 1.0e-15, ITER_max = 50000, ITER_min = 200, proposal_var = 0.01, technology = 0, H = 1000000000)
{

	MH_count_reject = 0
	Inner_count_reject = 0
	Z_list = matrix(0,ITER_max, k) 

	VAR_MH_old = 1000

	iter = 1
	GAP = 1000
		
	out_list = list()
	
	#---------------------------------------------
	# Feasible point:
	#---------------------------------------------

	A_kk = A_mat(1) 

	b_k = c(A_kk[1:m,kk], 0, 1)

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c_k[-kk] - c_kk

	c_k = c(max(c_k) - 0.9*c_k, 0)

	A_kk = cbind(A_kk[1:m,-kk] - A_kk[1:m,kk],  -A_kk[1:m,kk])
	
	AA_kk = rbind(A_kk, c(rep(1,k-1), -1), c(rep(0,k-1), 1))

	DIR = c(rep(">=", m), "<=", ">=")

	out_EDM = lp(direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_feasible = out_EDM$solution[1:(k-1)]
	delta_feasible = out_EDM$solution[k]
	
	#---------------------------------------------
	# Cost minimization:
	#---------------------------------------------

	c_k = c(P_I[kk,]%*%t(XX))

	c_k = c_k/max(abs(c_k))

	c_kk = c_k[kk]

	c_k = c(c_k[-kk] - c_kk, c_kk)

	out_EDM = lp (direction = "min", objective.in = c_k, const.mat = AA_kk, const.dir = DIR, const.rhs = b_k)
	z_opt = out_EDM$solution[1:(k-1)]
	delta_opt = out_EDM$solution[k]
		
	#################################################################
	#
	# Obtaining the expected posterior by MCMC
	#
	#################################################################

	set.seed(1234)

	# Sort activity vector in decreasing order
	
	if(max(z_opt) > 0){		
		H_set= order(z_opt, decreasing = TRUE)
	}else{		
		H_set= 1:(k-1)			
	}
	
	if(sum(abs(z_feasible - z_opt)) <= 0)
	{
		Z_list = matrix(rep(z_feasible, each = ITER_max), nrow = ITER_max, byrow = TRUE)
		
		z_plot = z_feasible[1]

		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)
		
		plot(0, pch = '')
		legend("topleft", legend = c(paste("Expectation = ", z_plot, sep = ''), paste("Optimal solution = ", z_opt[1], sep = '')), col=c(3, 2), lwd=3, cex=1)

		dev.off()
		
		out_list[[1]] = c(z_opt, delta_opt[1])		
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)	
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
	else
	{		
		delta_feasible = (0.99*delta_opt + 0.01*delta_feasible)
		z_opt = z_opt[H_set]
			
		# Initialization

		z_old = z_feasible
		delta_old = max(1, delta_feasible)
		
		sss = proposal_var

		while(((GAP >= tolerance) & (iter <= ITER_max)) || (iter < ITER_min) )
		{
			if((iter > 100) & ((iter/10) == round(iter/10)) )
			{			
				Rej_rate = MH_count_reject/(1 + iter)
				
				z_old = colMeans(Z_list[1:iter, (1:k-1)])
				delta_old = max(1, mean(Z_list[1:iter, k]))
												
				if(Rej_rate > 0.5)
				{ 
					sss = max(1.0e-10, sss*(1 - Rej_rate*0.01))
				}
				else
				{
					sss = max(1.0e-10, sss*(1 + (1-Rej_rate)*0.01))
				}
			}

			Simul = Simulate_z_NDRS(z_old, z_opt, delta_old, delta_opt, sss, A_kk, H_set, A_kk[1:m,kk])
			z_next = Simul[[1]]
			delta_next = Simul[[2]]
			Count_reject = Simul[[3]]
			
			Inner_count_reject = Inner_count_reject +  Count_reject

			if(Count_reject >= 1000)
			{
				prob_candidate_next = 1
				prob_candidate_old = 1
			}
			else
			{
				prob_candidate_next = g_delta_NDRS(z_next, z_old, delta_next, delta_old, sss) 
				prob_candidate_old = g_delta_NDRS(z_old, z_next, delta_old, delta_next, sss) 
			}
					
			AA = -H*as.numeric(c_k[H_set]%*%z_next + c_kk*delta_next) + H*as.numeric(c_k[H_set]%*%z_old + c_kk*delta_old) + prob_candidate_old - prob_candidate_next
			
			uu = runif(1)

			if(log(uu) <= AA)
			{
				z_old = z_next
				delta_old = delta_next
			}
			else
			{
				MH_count_reject = MH_count_reject + 1
			}

			Z_list[iter, ] = c(z_old, delta_old)
			
			if((iter > MH_count_reject + 10) & ((iter/500) == round(iter/500)))
			{
				VAR_MH = sqrt(diag(cov(Z_list[1:iter, ])))

				GAP = mean(abs(VAR_MH - VAR_MH_old))/(1 + max(VAR_MH_old))

				VAR_MH_old = VAR_MH

			}

			cat(c(round(iter), round(Inner_count_reject/(1 + iter), 3),  round(MH_count_reject/(1 + iter), 3), round(mean(Z_list[1:iter, 1]),3) , max(VAR_MH_old), GAP, '\n'), sep = "\t") 

			iter = iter + 1

		}

		ITER_max = iter -1
		
		BURN_IN = round(ITER_max/20)
		BURN_IN = 1
	
		filenamePLOT = paste("Plot_posterior_", kk, ".pdf", sep = "")

		pdf(filenamePLOT)

		MCMC_ind = seq(BURN_IN, ITER_max, by = 1)
		z_plot = Z_list[MCMC_ind, 1]
		Z_save = Z_list[MCMC_ind, ]
		
		z_down = min(z_plot)
		z_up = max(z_plot)
			
		hist(z_plot, breaks = 100, ylab = "Probability", xlab = "z[1]", xlim = c(z_down, z_up), main = "", axes = FALSE)
		legend("topleft", legend=c(paste("Expectation = ", mean(z_plot), sep = ''), paste("Optimal solution = ", z_opt[1], sep = ''), paste("Rejection rate = ", round(MH_count_reject/(1 + iter), 3), sep = '') ), col=c(3, 2, 0), lwd=3, cex=1)
		
		dev.off()		
		
		out_list[[1]] = Z_save
		out_list[[2]] = round(Inner_count_reject/(1 + iter), 3)
		out_list[[3]] = round(MH_count_reject/(1 + iter), 3)
		
		filenameSAVE = paste("Results_posterior_", kk, ".RData", sep = "")
		save(out_list, file = filenameSAVE)
		
	}
			
	return(out_list)


}



