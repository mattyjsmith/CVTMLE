******************************************************************************
* Title: Performance of Cross-Validated Targeted Maximum Likelihood Estimation 
* Author: Matthew J. Smith
******************************************************************************

*****	
* Change working directory
*****
	cd "YOUR WORKING DIRECTORY"
	clear

*****
* Preliminaries
*****
	set autotabgraphs on
	*set obs 1000
	*set seed 1
	local obs 1000
	local reps 1000

	
*****
* Run the Simulations
*****
	
/* 

Note! 

You will need to change the values for:
	1. The intercept (alpha) for the exposure model (to change the prevalence 
		of the exposure).
	2. The coefficient for the interaction (btreat2) in the outcome model (to 
		alter the extrapolation issue).

*/
	
	
	
* Create file to hold estimates from the simulation
	capture postclose ests
	postfile ests int(repno) float(ATE SE Method Delta p1theo p0theo p1 p0 logORtheo logOR SElogOR) using simests, replace
	
	* Start time of simulations
	scalar t1 = c(current_time)
	
	* Run the repetitions
	timer on 1
	qui {
		
		noi _dots 0, title("Simulation running...")
			
		forval rep = 1/`reps' {
			
			* Preliminaries
			clear
			set seed `rep' 
			set obs `obs'
			
			* Generate covariates
			gen z1 = rbinomial(1,0.1)
			gen z2 = rbinomial(1,0.4)
			gen z3 = rnormal(0,1)
			gen z4 = rbinomial(1,0.7)
			gen z5 = rbinomial(1,0.5)
			gen z6 = rnormal(0,1)
			gen z7 = rbinomial(1,0.3)
			gen z8 = rbinomial(1,0.8)
			gen z9 = rnormal(0,1)
			
			
			* Set parameters
			local alpha = -0.45  // Binary: 50% = -0.45, 80% = 1.05, Continuous: 50% = -0.35, 80% = 1.75
			local bz1 = log(5)     // Set 1.5 if binary has small effect, set 5 if large effect
			local bz6 = log(1.5)   // Set 1.5 if continuous has small effect, set 2.5 if continuous has large effect
			local btreat = log(1.75)	   // Coef for exposure
			local btreat2 = 0.0*`btreat'   // Coef for extrapolation issue
			
			* Exposure
			gen probA = invlogit(`alpha' + `bz1'*(z1) + log(1.5)*(z2) - log(1.5)*(z4) - `bz6'*(z6) + log(1.5)*(z7) + log(1.5)*(z8))			
			gen A = rbinomial(1, invlogit(`alpha' + `bz1'*(z1) + log(1.5)*(z2) - log(1.5)*(z4) - `bz6'*(z6) + log(1.5)*(z7) + log(1.5)*(z8)))
			*tab A
			
			* Outcome
				* btreat: log(1.75) or log(1) under H1 and H0, respectively. Coef for exposure
				* coeff.extrapol: 0.0, 0.3, 0.9, or 2.0. Coef for interaction between A and Z (outcome status)
				* t.theo: prob for this variable is either 0.8 or 0.5.	
			gen Y = rbinomial(1, invlogit(-0.8 + `btreat'*A + `btreat2'*A*z1 + log(1.5)*z1 + log(1.5)*z2 - log(1.5)*z3 - log(1.5)*z4 + log(1.5)*z5 + log(1.5)*z6))
			
			
			* Theo outcome
			sum probA
			local meanA r(mean) 
			gen ttheo = rbinomial(1, `meanA')
			gen ytheo = rbinomial(1, invlogit(-0.8 + `btreat'*ttheo + `btreat2'*ttheo*z1 + log(1.5)*z1 + log(1.5)*z2 - log(1.5)*z3 - log(1.5)*z4 + log(1.5)*z5 + log(1.5)*z6))
			
			
			// Run "True ATE"
			glm ytheo ttheo, f(b) link(logit)
			local logORtheo = e(b)[1,1]
			local p1theo = invlogit(e(b)[1,2] + e(b)[1,1])	// Logistic cumulative distribution function
			local p0theo = invlogit(e(b)[1,2])				// Logistic cumulative distribution function
			local Delta = `p1theo' - `p0theo'
			di `Delta'
			
			
			// Run TMLE
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, tmle elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 1		
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')
			*/
				
			
			
			// Run CVTMLE
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, cvtmle cvfolds(10) elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 2
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')
			
			
			// Run CVTMLE(Qg)
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, cvtmleQg cvfolds(10) elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 3
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')
			
			
			// Run TMLE with RF
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, tmleglsrf elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 4	
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')
		
			
			
			// Run CVTMLE(Q) with RF
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, cvtmleglsrf cvfolds(10) elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 5	
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')

			
			
			
			// Run CVTMLE(Qg) with RF
			*use "simCVeltmledata", clear
			eltmle Y A z1 z2 z3 z4 z5 z6, cvtmleQgglsrf cvfolds(10) elements
			local logOR = log(r(MOR))
			local SElogOR = r(SE_log_MOR)
			local ATE = r(ATEtmle)
			local SE = r(ATE_SE_tmle)
			local Method = 6
			sum _POM1
			local p1 = r(mean)
			sum _POM0
			local p0 = r(mean)
				* Store the estimates
				post ests (`rep') (`ATE') (`SE') (`Method') (`Delta') (`p1theo') (`p0theo') (`p1') (`p0') (`logORtheo') (`logOR') (`SElogOR')
			
			
		
		* Dot for completion of rep
		noi _dots `rep' 0
		
		}
	}
	timer off 1
	postclose ests
	
	* End time of simulations
	scalar t2 = c(current_time)
	
	* Computational time
	display (clock(t2, "hms") - clock(t1, "hms")) / 1000 " seconds"	
	display ((clock(t2, "hms") - clock(t1, "hms")) / 1000)/(60*60) " hours"	
	
		
	* load estimates data
		use "simests.dta", clear
	
	
	* Performance measures
		* ATE
		qui: sum Delta
		local Deltaest = r(mean)
		qui: sum ATE
		local trueATE = r(mean)
		di "Relative bias = " (abs(`Deltaest' - `trueATE')/`trueATE')*100 "%"
		simsum ATE, true(`Deltaest') meth(Method) id(repno) se(SE)
		

		
		
		
*****
* Produce graphs
*****


	cd "YOUR WORKING DIRECTORY"
	clear
	
	*set autotabgraphs on
	set obs 1000
	set seed 1
	
	* Generate covariates
	gen z1 = rbinomial(1,0.1)
	gen z2 = rbinomial(1,0.4)
	gen z3 = rnormal(0,1)
	gen z4 = rbinomial(1,0.7)
	gen z5 = rbinomial(1,0.5)
	gen z6 = rnormal(0,1)			// z6 is continuous variable
	gen z7 = rbinomial(1,0.3)
	gen z8 = rbinomial(1,0.8)
	gen z9 = rnormal(0,1)
			
			
	* Specify parameters
	local alpha = -0.45			// -0.45 = 50%, 1.05 = 80%
	local bz = log(5)			// Coef for binary variable
	local bz6 = log(1.5)		// Coef for positivity covariate in exposure model
	local btreat = log(1.75) 	// Coef for exposure in outcome model.  log(1.75) or log(1) under H1 and H0, respectively.
	local btreat2 = 0.0*`btreat' // Coef for extrapolation issue. 0.0, 0.3, 0.9, or 2.0. Coef for interaction between A and Z (outcome status).
		
		
	* Generate exposure
	capture drop probA* A* 
	
	gen probA = invlogit(`alpha' + `bz'*(z1) + log(1.5)*(z2) - log(1.5)*(z4) - ///
							`bz6'*(z6) + log(1.5)*(z7) + log(1.5)*(z8))
	gen A = rbinomial(1, invlogit(`alpha' + `bz'*(z1) + log(1.5)*(z2) - log(1.5)*(z4) - ///
							`bz6'*(z6) + log(1.5)*(z7) + log(1.5)*(z8)))
	tab A
	
		
	* Observed outcome			
	capture drop Y* probY*	
	
	gen probY = invlogit(-0.8 + `btreat'*A + `btreat2'*A*z1 + log(1.5)*z1 + ///
							log(1.5)*z2 - log(1.5)*z3 - log(1.5)*z4 + log(1.5)*z5 + log(1.5)*z6)
	gen Y = rbinomial(1, invlogit(-0.8 + `btreat'*A + `btreat2'*A*z1 + log(1.5)*z1 + ///
							log(1.5)*z2 - log(1.5)*z3 - log(1.5)*z4 + log(1.5)*z5 + log(1.5)*z6))
	tab Y
	
		
	* Theo outcome
	capture drop ttheo* ytheo*		
	sum probA
	local meanA r(mean) 
	gen ttheo = rbinomial(1, `meanA')
	gen ytheo = rbinomial(1, invlogit(-0.8 + `btreat'*ttheo + `btreat2'*ttheo*z1 + log(1.5)*z1 + log(1.5)*z2 - log(1.5)*z3 - log(1.5)*z4 + log(1.5)*z5 + log(1.5)*z6))
	

	
	
*****
* Overlap plots
*****

	* Overlap
		teffects ipw (Y) (A z1 z2 z3 z4 z5 z6)
		teoverlap, title("P[A = 1] = 0.5") ///
			 legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0)) ///
			 yscale(r(0 8)) ylabel(0(1)8) xscale(r(0 1)) xlabel(0(0.2)1) ///
			 name(Overlap50, replace) saving(Overlap50, replace)
		

	* Graph the outcome 
		twoway	(lowess probY z1 if A==1, lcolor(red)  lpattern(solid) xlabel(0(1)1) ///
												yscale(r(0 1)) ylabel(0(0.1)1)) ///
				(lowess probY z1 if A==0, lcolor(blue) lpattern(solid) xlabel(0(1)1) ///
												yscale(r(0 1)) ylabel(0(0.1)1)), ///
					legend(order(1 "Treated" 2 "Not treated") position(11) cols(1) ring(0)) ///	
					ytitle("P(Y=1 | A,z1)") ///
					title("(A)") ///
					name(PYNone50, replace)
					
		
		
					
		

		
		
		
		