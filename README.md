Supplementary Materials to Swihart & Bandyopadhyay (2020)
================

**Note: the dataset in this repo is a subset of the one used in the
paper. Please contact Dipankar Bandyopadhyay (<Bandyop@vcuhealth.org>)
for full data**

## Conditional Proportional Hazards (Recursive-*Ω* in R)

``` r
start_time <- Sys.time()
cphz_recomega <- parfm(Surv(time, status) ~
                               diabetes_yn_10 +
                               tobacco_use_yn_10 +
                               male_yn_10 +
                               age_stdzd +
                               bop_stdzd +
                               plaque_stdzd +
                               calmean_stdzd +
                               crown_yn_10 +
                               molar_yn_10,
                             cluster="id", 
                             data=creighton_sub, 
                             dist="weibull", 
                             frailty="possta")
end_time <- Sys.time()
cphz_recomega_time <- end_time - start_time
## put parfm results in 95% CI format
parfm_results <-
data.frame(value=cphz_recomega[,1], 
           low95=cphz_recomega[,1]-1.96*cphz_recomega[,2], 
           upp95=cphz_recomega[,1]+1.96*cphz_recomega[,2],
           loglik=sprintf("%.10f",attr(cphz_recomega, "loglik"))
           )
```

## Conditional Proportional Hazards (Stirling-Static in R)

``` r
start_time <- Sys.time()
## set up integrated likelihood:
integrated_likelihood <- function(parms){
  
        a  <- parms[ 1];
  shape_in <- parms[ 2];
  beta00   <- parms[ 3];
  beta01   <- parms[ 4];
  beta02   <- parms[ 5];
  beta03   <- parms[ 6];
  beta04   <- parms[ 7];
  beta05   <- parms[ 8];
  beta06   <- parms[ 9];
  beta07   <- parms[10];
  beta08   <- parms[11];
  beta09   <- parms[12];
  
nll <- 
  -sum(
    log(
      creighton_sub[,
           {
             eta <- beta00 + 
               beta01*diabetes_yn_10 +
               beta02*tobacco_use_yn_10 +
               beta03*male_yn_10 +
               beta04*age_stdzd +
               beta05*bop_stdzd +
               beta06*plaque_stdzd +
               beta07*calmean_stdzd +
               beta08*crown_yn_10 +
               beta09*molar_yn_10;          
                          
               shape <- shape_in;            ## COND _PH shape     
               scale <- exp(eta);            ## COND _PH scale
             
            #  shape <- shape_in         ;   ## COND AFT shape
            #  scale <- exp(-shape * eta);   ## COND AFT scale
              
            #  shape <- shape_in/a;          ## MARG _PH shape
            #  scale <- exp(eta/a);          ## MARG _PH scale
              
            #  shape <- shape_in/a         ; ## MARG AFT shape     
            #  scale <- exp(-shape * eta/a); ## MARG AFT scale

               prod(scale^status) *
               shape^sum(status) *
               prod(time^status)^(shape-1) *
               integral_term(sum_status = sum(status), a, s=sum(scale*time^shape) )
               

           },
           by=id]$V1
   )
)

nll
}

start_parms_val <- c(0.5, 1, 0.10, rep(0.1, 9))

cphz_stirling <-
optim(start_parms_val,
      integrated_likelihood,
      method="L-BFGS-B",
      
      lower=c(0.1,0.1, -7, rep(-3, 9)),  ## COND  PH lower bounds     
      upper=c(0.9,2.0,  3, rep( 3, 9)),  ## COND  PH upper bounds     

      hessian=TRUE,
      control=list(factr=10^3.9))

#sprintf("%.10f", -cphz_stirling$value)
#sprintf("%.10f",attr(cphz_recomega, "loglik"))

## YES!
fit<-cphz_stirling
fisher_info <- solve(fit$hessian)
prop_sigma<-matrix(sqrt(diag(fisher_info)))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
cphz_stirling_results <-data.frame(value=fit$par, low95=lower, upp95=upper, loglik=sprintf("%.10f", -cphz_stirling$value))

row.names(cphz_stirling_results) <- row.names(parfm_results)
rownames(cphz_stirling_results)[1] <- "alpha = (1-nu)"
rownames(cphz_stirling_results)[3] <- "log(lambda)"

end_time <- Sys.time()
cphz_stirling_time <- end_time - start_time
```

## Compare the two closed forms:

-   Conditional Proportional Hazards - Recursive-*Ω* time: 12.89 seconds
-   Conditional Proportional Hazards - Static-Stirling time: 3.32
    seconds

Static-Stirling is faster. Also provided the same likelihood:

-   Conditional Proportional Hazards - Recursive-*Ω* likelihood:
    -100.5676297313
-   Conditional Proportional Hazards - Static-Stirling likelihood:
    -100.5676297294

The estimates and 95%CI are very similar as well:

-   Conditional Proportional Hazards - Recursive-*Ω* estimates:

``` r
  print(parfm_results[-1*c(1,2,3),-4], digits = 3)
```

    ##                     value  low95 upp95
    ## diabetes_yn_10     1.6826 -0.218 3.583
    ## tobacco_use_yn_10  0.1168 -1.370 1.603
    ## male_yn_10        -0.3110 -2.006 1.384
    ## age_stdzd         -0.9324 -2.019 0.154
    ## bop_stdzd         -0.0743 -0.459 0.310
    ## plaque_stdzd      -0.3815 -1.116 0.353
    ## calmean_stdzd      0.3230 -0.149 0.795
    ## crown_yn_10        2.0242  1.040 3.008
    ## molar_yn_10        0.0492 -0.790 0.889

-   Conditional Proportional Hazards - Static-Stirling estimates:

``` r
  print(cphz_stirling_results[-1*c(1,2,3),-4], digits=3)  
```

    ##                     value  low95 upp95
    ## diabetes_yn_10     1.6826 -0.221 3.586
    ## tobacco_use_yn_10  0.1168 -1.397 1.630
    ## male_yn_10        -0.3110 -2.020 1.398
    ## age_stdzd         -0.9324 -2.019 0.154
    ## bop_stdzd         -0.0743 -0.460 0.312
    ## plaque_stdzd      -0.3815 -1.116 0.353
    ## calmean_stdzd      0.3230 -0.150 0.796
    ## crown_yn_10        2.0242  1.029 3.019
    ## molar_yn_10        0.0492 -0.793 0.891

You’ll note that the first three rows are omitted. They are the same,
but displayed in different ways. `parfm` is parameterized as `nu` which
equals Kendall’s Tau which is 1-alpha in the Static Stirling
formulation. Also, `parfm` gives lambda, whereas Static-Stirling give
log(lambda):

-   Conditional Proportional Hazards - Recursive-*Ω* estimates:

``` r
  print(parfm_results[1*c(1,2,3),-4], digits = 3)
```

    ##          value    low95  upp95
    ## nu     0.28082  0.00879 0.5529
    ## rho    0.77964  0.52919 1.0301
    ## lambda 0.00981 -0.00615 0.0258

-   Conditional Proportional Hazards - Static-Stirling estimates:

``` r
  print(cphz_stirling_results[1*c(1,2,3),-4], digits=3)  
```

    ##                 value  low95  upp95
    ## alpha = (1-nu)  0.719  0.446  0.992
    ## rho             0.780  0.529  1.031
    ## log(lambda)    -4.624 -6.360 -2.888

Sure, some of the estimates are slightly off, and this discrepancy
diminishes as the sample size grows. It is best to focus on how similar
the likelihood estimates are.

At this point, we have estimated the conditional proportional hazards
perspective-parameterization. The conditional proportional hazards
perspective-parameterization is the sole perspective-parameterization
offered currently by `parfm`. Exponentiating the coefficients give the
cluster-specific proportional hazard ratios and exponentiating the 95%
CI bounds give the 95% CI bounds for the cluster-specific hazard ratios.

But what if we wanted the other 3 perspective-parameterizations? If we
only wanted estimates for marginal proportional hazards, conditional
acceleration factor, and/or marginal acceleration factor, we just need
to take the direct esimates of the conditional proportional hazards
perspective-parameterization and do some calculations as per this table:

![see Table 3 in the
paper](table3.png "Note: The diagonal has the direct estimates for each perspective-parameterizations. This is Table 3 from the paper.")

Our estimate of *α* is 0.719, the estimate of
*β*<sub>*h*</sub><sup>*c*</sup> is 1.683 for diabetes, and the estimate
of the shape *p* is 0.780.

So from the 1st row of Table 3 above we can cal

-   Conditional hazard ratios: exp (*β*<sub>*h*</sub><sup>*c*</sup>) =
    exp (1.683) = 5.379

-   Marginal hazard ratios: exp (*α* × *β*<sub>*h*</sub><sup>*c*</sup>)
    = exp (0.719 × 1.683) = 3.354

-   Conditional acceleration factor:
    $\\exp(\\frac{1}{-p} \\times \\beta^c\_h)$ = exp (1.683) = 5.379

-   Marginal acceleration factor:
    $\\exp(\\frac{\\alpha}{-p} \\times \\beta^c\_h)$ =
    exp (0.719 × 1.683) = 3.354
