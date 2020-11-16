Supplementary Materials to Swihart & Bandyopadhyay (2020)
================

**Note: the dataset in this repo is a subset of the one used in the
paper. Please contact Dipankar Bandyopadhyay (<Bandyop@vcuhealth.org>)
for full data**

``` r
knitr::opts_chunk$set(echo = TRUE)
options(digits=3)
source("./R/integral_terms_all_v02.R")
creighton_sub <- readRDS("creighton_sub_v01.RDS")
library(broman)
library(data.table)
library(parfm)
```

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

-   Conditional Proportional Hazards - Recursive-*Ω* time: 12.47 seconds
-   Conditional Proportional Hazards - Static-Stirling time: 3.37
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

So from the 1st row of Table 3 above we can calculate from the direct
estimates all the other estimands of interest:

-   Conditional hazard ratios: exp (*β*<sub>*h*</sub><sup>*c*</sup>) =
    exp (1.683) = 5.379

-   Marginal hazard ratios: exp (*α* × *β*<sub>*h*</sub><sup>*c*</sup>)
    = exp (0.719 × 1.683) = 3.354

-   Conditional acceleration factor:
    exp (*β*<sub>*h*</sub><sup>*c*</sup>/( − *p*)) =
    exp (1.683/( − 0.780)) = 0.116

-   Marginal acceleration factor:
    exp ((*α* × *β*<sub>*h*</sub><sup>*c*</sup>)/( − *p*)) =
    exp ((0.719 × 1.683)/( − 0.780)) = 0.212

If we want the uncertainty of the other esitmands of interest (i.e., the
standard error of the Marginal acceleration factor or the 95% CIs of the
Marginal hazard ratios, etc) then we need to do one of two things:

1.  Derive some sort of Delta-Method type variance for the estimand of
    interest that wasn’t directly estimated
2.  Directly estimate the estimand of interest! That is, rerun the same
    model but change the `shape` and `scale` in the negative
    log-likelihood (`nll`) in the code above.

We choose option 2! It’s easy. Let’s do it for Marginal Hazard Ratios.

## Marginal Proportional Hazards (Stirling-Static in R)

We treat the conditional proportional hazards as a “seed” model – that
is, run it first. It is the most forgiving in term of initial starting
values and boundaries. Once estimates are obtained from it, fitting the
other perspective-parameterizations can use adjusted starting values and
boundaries. For instance, we know that our estimate of *α* is 0.719 and
that won’t change for any of the other perspetive-parameterization
combinations. We also know from the theory that all the betas will be
multipled by alpha, so we can adjust the bounds by multiplying by alpha.

``` r
start_time <- Sys.time()
## set up integrated likelihood: uncomment the shape, scale you want, comment the rest
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
                          
            #  shape <- shape_in;            ## COND _PH shape     
            #  scale <- exp(eta);            ## COND _PH scale
             
            #  shape <- shape_in         ;   ## COND AFT shape
            #  scale <- exp(-shape * eta);   ## COND AFT scale
              
               shape <- shape_in/a;          ## MARG _PH shape
               scale <- exp(eta/a);          ## MARG _PH scale
              
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

start_parms_val <- c(cphz_stirling$par[1],
                     cphz_stirling$par[c(-1)]*cphz_stirling$par[1])*1.05

mphz_stirling <-
optim(start_parms_val,
      integrated_likelihood,
      method="L-BFGS-B",
      
      lower=start_parms_val*(0.80*(start_parms_val>0) + 1.20*(start_parms_val<0)),  ## MARG  PH lower bounds     
      upper=start_parms_val*(1.20*(start_parms_val>0) + 0.80*(start_parms_val<0)),  ## MARG  PH upper bounds     

      hessian=TRUE,
      control=list(factr=10^3.9))

## YES!
fit<-mphz_stirling
fisher_info <- solve(fit$hessian)
prop_sigma<-matrix(sqrt(diag(fisher_info)))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
mphz_stirling_results <-data.frame(value=fit$par, low95=lower, upp95=upper, loglik=sprintf("%.10f", -mphz_stirling$value))

row.names(mphz_stirling_results) <- row.names(parfm_results)
rownames(mphz_stirling_results)[1] <- "alpha = (1-nu)"
rownames(mphz_stirling_results)[3] <- "log(lambda)"

end_time <- Sys.time()
mphz_stirling_time <- end_time - start_time
```

## Compare the two closed forms:

-   Conditional Proportional Hazards - Recursive-*Ω* time: 12.47 seconds
-   Marginal Proportional Hazards - Static-Stirling time: 2.53 seconds

Static-Stirling is faster. Also provided the same likelihood:

-   Conditional Proportional Hazards - Recursive-*Ω* likelihood:
    -100.5676297313
-   Marginal Proportional Hazards - Static-Stirling likelihood:
    -100.5676297295

The estimates and 95%CI are very similar as well – just multiply *α* to
get the marginal:

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

-   Marginal Proportional Hazards - Static-Stirling estimates:

``` r
  print(mphz_stirling_results[-1*c(1,2,3),-4], digits=3)  
```

    ##                     value  low95 upp95
    ## diabetes_yn_10     1.2101 -0.259 2.679
    ## tobacco_use_yn_10  0.0840 -1.002 1.170
    ## male_yn_10        -0.2236 -1.464 1.017
    ## age_stdzd         -0.6706 -1.518 0.177
    ## bop_stdzd         -0.0534 -0.329 0.222
    ## plaque_stdzd      -0.2744 -0.789 0.240
    ## calmean_stdzd      0.2323 -0.118 0.582
    ## crown_yn_10        1.4558  0.598 2.313
    ## molar_yn_10        0.0354 -0.569 0.640

You’ll note that the first three rows are omitted. *α* remains the same,
but shape changes.

-   Conditional Proportional Hazards - Recursive-*Ω* estimates:

``` r
  print(parfm_results[1*c(1,2,3),-4], digits = 3)
```

    ##          value    low95  upp95
    ## nu     0.28082  0.00879 0.5529
    ## rho    0.77964  0.52919 1.0301
    ## lambda 0.00981 -0.00615 0.0258

-   Marginal Proportional Hazards - Static-Stirling estimates:

``` r
  print(mphz_stirling_results[1*c(1,2,3),-4], digits=3)  
```

    ##                 value  low95  upp95
    ## alpha = (1-nu)  0.719  0.446  0.992
    ## rho             0.561  0.302  0.820
    ## log(lambda)    -3.326 -4.769 -1.882

## Conditional Accerlation Factor (Stirling-Static in R)

We treat the conditional proportional hazards as a “seed” model – that
is, run it first. It is the most forgiving in term of initial starting
values and boundaries. Once estimates are obtained from it, fitting the
other perspective-parameterizations can use adjusted starting values and
boundaries. For instance, we know that our estimate of *α* is 0.719 and
that won’t change for any of the other perspetive-parameterization
combinations. We also know from the theory that all the betas will be
multiplied by -1/p, so we can adjust the bounds and starting values
accordingly.

``` r
start_time <- Sys.time()
## set up integrated likelihood: uncomment the shape, scale you want, comment the rest
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
                          
            #  shape <- shape_in;            ## COND _PH shape     
            #  scale <- exp(eta);            ## COND _PH scale
             
               shape <- shape_in         ;   ## COND AFT shape
               scale <- exp(-shape * eta);   ## COND AFT scale
              
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

start_parms_val <- c( cphz_stirling$par[1],
                      cphz_stirling$par[2],
                     -cphz_stirling$par[c(-1,-2)]/cphz_stirling$par[2])

caft_stirling <-
optim(start_parms_val,
      integrated_likelihood,
      method="L-BFGS-B",
      
      lower=start_parms_val*(0.80*(start_parms_val>0) + 1.20*(start_parms_val<0)),  ## MARG  PH lower bounds     
      upper=start_parms_val*(1.20*(start_parms_val>0) + 0.80*(start_parms_val<0)),  ## MARG  PH upper bounds     

      hessian=TRUE,
      control=list(factr=10^3.9))


## YES!
fit<-caft_stirling
fisher_info <- solve(fit$hessian)
prop_sigma<-matrix(sqrt(diag(fisher_info)))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
caft_stirling_results <-data.frame(value=fit$par, low95=lower, upp95=upper, loglik=sprintf("%.10f", -caft_stirling$value))

row.names(caft_stirling_results) <- row.names(parfm_results)
rownames(caft_stirling_results)[1] <- "alpha = (1-nu)"
rownames(caft_stirling_results)[3] <- "log(lambda)"

end_time <- Sys.time()
caft_stirling_time <- end_time - start_time
```

## Compare the two closed forms:

-   Conditional Proportional Hazards - Recursive-*Ω* time: 12.47 seconds
-   Conditional Acceleration Factor - Static-Stirling time: 0.97 seconds

Static-Stirling is faster. Also provided the same likelihood:

-   Conditional Proportional Hazards - Recursive-*Ω* likelihood:
    -100.5676297313
-   Conditional Acceleration Factor- Static-Stirling likelihood:
    -100.5676297294

The estimates and 95%CI are very similar as well – just divide  − *p* to
get them:

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

-   Conditional Acceleration Factor - Static-Stirling estimates:

``` r
  print(caft_stirling_results[-1*c(1,2,3),-4], digits=3)  
```

    ##                     value  low95  upp95
    ## diabetes_yn_10    -2.1582 -4.631  0.315
    ## tobacco_use_yn_10 -0.1498 -2.089  1.790
    ## male_yn_10         0.3988 -1.800  2.597
    ## age_stdzd          1.1960 -0.172  2.564
    ## bop_stdzd          0.0953 -0.398  0.589
    ## plaque_stdzd       0.4894 -0.472  1.451
    ## calmean_stdzd     -0.4143 -1.021  0.192
    ## crown_yn_10       -2.5963 -3.999 -1.194
    ## molar_yn_10       -0.0631 -1.145  1.019

You’ll note that the first three rows are omitted. *α* and the shape
remains the same. the shape only changes between marginal and
conditional perspectives.

-   Conditional Proportional Hazards - Recursive-*Ω* estimates:

``` r
  print(parfm_results[1*c(1,2,3),-4], digits = 3)
```

    ##          value    low95  upp95
    ## nu     0.28082  0.00879 0.5529
    ## rho    0.77964  0.52919 1.0301
    ## lambda 0.00981 -0.00615 0.0258

-   Conditional Acceleration Factor - Static-Stirling estimates:

``` r
  print(caft_stirling_results[1*c(1,2,3),-4], digits=3)  
```

    ##                value low95 upp95
    ## alpha = (1-nu) 0.719 0.446 0.992
    ## rho            0.780 0.529 1.031
    ## log(lambda)    5.931 3.329 8.534

## Marginal Accerlation Factor (Stirling-Static in R)

We treat the conditional proportional hazards as a “seed” model – that
is, run it first. It is the most forgiving in term of initial starting
values and boundaries. Once estimates are obtained from it, fitting the
other perspective-parameterizations can use adjusted starting values and
boundaries. For instance, we know that our estimate of *α* is 0.719 and
that won’t change for any of the other perspetive-parameterization
combinations. We also know from the theory that all the betas will be
multiplied by -*α*/p, so we can adjust the bounds and starting values
accordingly.

``` r
start_time <- Sys.time()
## set up integrated likelihood: uncomment the shape, scale you want, comment the rest
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
                          
            #  shape <- shape_in;            ## COND _PH shape     
            #  scale <- exp(eta);            ## COND _PH scale
             
            #  shape <- shape_in         ;   ## COND AFT shape
            #  scale <- exp(-shape * eta);   ## COND AFT scale
              
            #  shape <- shape_in/a;          ## MARG _PH shape
            #  scale <- exp(eta/a);          ## MARG _PH scale
              
               shape <- shape_in/a         ; ## MARG AFT shape     
               scale <- exp(-shape * eta/a); ## MARG AFT scale

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

start_parms_val <- c( cphz_stirling$par[1],
                      c(cphz_stirling$par[2],
                        -cphz_stirling$par[c(-1,-2)]/cphz_stirling$par[2]) * cphz_stirling$par[1]
                      )


maft_stirling <-
optim(start_parms_val,
      integrated_likelihood,
      method="L-BFGS-B",
      
      lower=start_parms_val*(0.80*(start_parms_val>0) + 1.20*(start_parms_val<0)),  ## MARG  PH lower bounds     
      upper=start_parms_val*(1.20*(start_parms_val>0) + 0.80*(start_parms_val<0)),  ## MARG  PH upper bounds     

      hessian=TRUE,
      control=list(factr=10^3.9))


## YES!
fit<-maft_stirling
fisher_info <- solve(fit$hessian)
prop_sigma<-matrix(sqrt(diag(fisher_info)))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
maft_stirling_results <-data.frame(value=fit$par, low95=lower, upp95=upper, loglik=sprintf("%.10f", -maft_stirling$value))

row.names(maft_stirling_results) <- row.names(parfm_results)
rownames(maft_stirling_results)[1] <- "alpha = (1-nu)"
rownames(maft_stirling_results)[3] <- "log(lambda)"

end_time <- Sys.time()
maft_stirling_time <- end_time - start_time
```

## Compare the two closed forms:

-   Conditional Proportional Hazards - Recursive-*Ω* time: 12.47 seconds
-   Marginal Acceleration Factor - Static-Stirling time: 1.08 seconds

Static-Stirling is faster. Also provided the same likelihood:

-   Conditional Proportional Hazards - Recursive-*Ω* likelihood:
    -100.5676297313
-   Marginal Acceleration Factor- Static-Stirling likelihood:
    -100.5676297294

The estimates and 95%CI are very similar as well – just multiply
 − *α*/*p* to get them:

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

-   Marginal Acceleration Factor - Static-Stirling estimates:

``` r
  print(maft_stirling_results[-1*c(1,2,3),-4], digits=3)  
```

    ##                     value  low95  upp95
    ## diabetes_yn_10    -1.5521 -3.479  0.375
    ## tobacco_use_yn_10 -0.1077 -1.500  1.285
    ## male_yn_10         0.2868 -1.310  1.883
    ## age_stdzd          0.8601 -0.222  1.942
    ## bop_stdzd          0.0685 -0.284  0.421
    ## plaque_stdzd       0.3519 -0.325  1.029
    ## calmean_stdzd     -0.2980 -0.750  0.154
    ## crown_yn_10       -1.8672 -3.094 -0.640
    ## molar_yn_10       -0.0454 -0.822  0.731

You’ll note that the first three rows are omitted. *α* and the shape
remains the same. the shape only changes between marginal and
conditional perspectives.

-   Conditional Proportional Hazards - Recursive-*Ω* estimates:

``` r
  print(parfm_results[1*c(1,2,3),-4], digits = 3)
```

    ##          value    low95  upp95
    ## nu     0.28082  0.00879 0.5529
    ## rho    0.77964  0.52919 1.0301
    ## lambda 0.00981 -0.00615 0.0258

-   Marginal Acceleration Factor - Static-Stirling estimates:

``` r
  print(maft_stirling_results[1*c(1,2,3),-4], digits=3)  
```

    ##                value low95 upp95
    ## alpha = (1-nu) 0.719 0.446 0.992
    ## rho            0.561 0.302 0.820
    ## log(lambda)    4.266 2.030 6.501

## FAQs:

### Q: What if my application needs an integrated term for n &gt; 32?

### A: Run this code and store the output in a file similar to `integral_terms_all_v02.R`. Our application needed `n in c(1:32)`. Remember for n=0 the integrated term is `exp(-s^a)`. Just start with 1 below and go to the max number of observed failures in the a cluster for your dataset. For instance, a cluster-randomized clinical trial with clusters of a size 100 would be covered by `n in c(1:100)`. We show the code for `n in c(1:7)` and the output, which can be copy-and-pasted into functions.

    ## copy output into functions like we did for integral_terms_all_v02.R
    library(copula)
    for(n in c(1:7)){
      term_string_collector <- NULL
      
      for(i in 0:(n-1)){
        
        within_parenthetical_coef <- c(1,rep(c(-1,1), length=(i)))  * rev(copula::Stirling2.all(i+1))
        
        term_string <-
          noquote(

            paste0("(",
                   paste0("(", within_parenthetical_coef,"*sa^", (i:0), ")", collapse="+"),
                   ")"
            )

          )
        
        term_string_collector <- c(term_string, term_string_collector)
      }

      term_string_collector
      

      
      n_gt1_rep <- 
        noquote(
          paste0("(",
                 paste0(rev(abs(copula::Stirling1.all(n))) , "*a^",((n-1):0), "*",term_string_collector, collapse="+"),
                 ")"
          )
        )
      
      print(paste0("This is n=", n, ":"))
      print(noquote(paste0( paste0("a*esa*s^(a-",n,")*"), n_gt1_rep)))
      print(noquote("/end"))
    }

The output:

    [1] "This is n=1:"
    [1] a*esa*s^(a-1)*(1*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=2:"
    [1] a*esa*s^(a-2)*(1*a^1*((1*sa^1)+(-1*sa^0))+1*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=3:"
    [1] a*esa*s^(a-3)*(1*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+3*a^1*((1*sa^1)+(-1*sa^0))+2*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=4:"
    [1] a*esa*s^(a-4)*(1*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+6*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+11*a^1*((1*sa^1)+(-1*sa^0))+6*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=5:"
    [1] a*esa*s^(a-5)*(1*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+10*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+35*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+50*a^1*((1*sa^1)+(-1*sa^0))+24*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=6:"
    [1] a*esa*s^(a-6)*(1*a^5*((1*sa^5)+(-15*sa^4)+(65*sa^3)+(-90*sa^2)+(31*sa^1)+(-1*sa^0))+15*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+85*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+225*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+274*a^1*((1*sa^1)+(-1*sa^0))+120*a^0*((1*sa^0)))
    [1] /end
    [1] "This is n=7:"
    [1] a*esa*s^(a-7)*(1*a^6*((1*sa^6)+(-21*sa^5)+(140*sa^4)+(-350*sa^3)+(301*sa^2)+(-63*sa^1)+(1*sa^0))+21*a^5*((1*sa^5)+(-15*sa^4)+(65*sa^3)+(-90*sa^2)+(31*sa^1)+(-1*sa^0))+175*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+735*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+1624*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+1764*a^1*((1*sa^1)+(-1*sa^0))+720*a^0*((1*sa^0)))
    [1] /end

Which can then be made into functions:

    integral_term_0001 <- function(esa,sa,s,a){
      a*esa*s^(a-1)*(1*a^0*((1*sa^0)))
    } 

    integral_term_0002 <- function(esa,sa,s,a){
      a*esa*s^(a-2)*(1*a^1*((1*sa^1)+(-1*sa^0))+1*a^0*((1*sa^0)))
    } 

    integral_term_0003 <- function(esa,sa,s,a){
      a*esa*s^(a-3)*(1*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+3*a^1*((1*sa^1)+(-1*sa^0))+2*a^0*((1*sa^0)))
    } 

    integral_term_0004 <- function(esa,sa,s,a){
      a*esa*s^(a-4)*(1*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+6*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+11*a^1*((1*sa^1)+(-1*sa^0))+6*a^0*((1*sa^0)))  
    } 

    integral_term_0005 <- function(esa,sa,s,a){
      a*esa*s^(a-5)*(1*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+10*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+35*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+50*a^1*((1*sa^1)+(-1*sa^0))+24*a^0*((1*sa^0)))  
    } 

    integral_term_0006 <- function(esa,sa,s,a){
      a*esa*s^(a-6)*(1*a^5*((1*sa^5)+(-15*sa^4)+(65*sa^3)+(-90*sa^2)+(31*sa^1)+(-1*sa^0))+15*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+85*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+225*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+274*a^1*((1*sa^1)+(-1*sa^0))+120*a^0*((1*sa^0)))  
    } 

    integral_term_0007 <- function(esa,sa,s,a){
      a*esa*s^(a-7)*(1*a^6*((1*sa^6)+(-21*sa^5)+(140*sa^4)+(-350*sa^3)+(301*sa^2)+(-63*sa^1)+(1*sa^0))+21*a^5*((1*sa^5)+(-15*sa^4)+(65*sa^3)+(-90*sa^2)+(31*sa^1)+(-1*sa^0))+175*a^4*((1*sa^4)+(-10*sa^3)+(25*sa^2)+(-15*sa^1)+(1*sa^0))+735*a^3*((1*sa^3)+(-6*sa^2)+(7*sa^1)+(-1*sa^0))+1624*a^2*((1*sa^2)+(-3*sa^1)+(1*sa^0))+1764*a^1*((1*sa^1)+(-1*sa^0))+720*a^0*((1*sa^0)))  
    } 

Which can then be called by a parent function:

    ## integral_term() based on a switch function.
    ## for each `n`, an integral term function

    integral_term <- function(sum_status,a,s){
      
      n <- sum_status
      
      sa <- s^a
      esa <- exp(-sa)
      terms <- NULL
      
      offset_1 <- sum_status+1
      switch(offset_1,
             
             exp(-s^a),
             
             integral_term_0001(esa,sa,s,a),
             integral_term_0002(esa,sa,s,a),
             integral_term_0003(esa,sa,s,a),
             integral_term_0004(esa,sa,s,a),
             integral_term_0005(esa,sa,s,a),
             integral_term_0006(esa,sa,s,a),
             integral_term_0007(esa,sa,s,a)

      )
    }



### Q: What about SAS?

### A: SAS! Computation times will take longer. Examples of the code are provided in the SAS folder in this repo. If at all possible, we’d recommend fitting the conditional proportional hazards model in R (with say, `parfm`) to get an idea of the estimates and bounds first. All the options used in the SAS code were deemed necessary for convergence/speed.
