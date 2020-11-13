---
title: "Supplementary Materials to Swihart & Bandyopadhyay (2020)"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---



**Note:  the dataset in this repo is a subset of the one used in the paper.  Please contact Dipankar Bandyopadhyay (Bandyop@vcuhealth.org) for full data**


## Conditional Proportional Hazards (Recursive-$\Omega$ in R)


```r
start_time <- Sys.time()
start_time
```

```
## [1] "2020-11-13 15:40:07 EST"
```

```r
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
end_time - start_time
```

```
## Time difference of 13.66609 secs
```

```r
cphz_recomega
```

```
## 
## Frailty distribution: positive stable 
## Baseline hazard distribution: Weibull 
## Loglikelihood: -100.568 
## 
##                   ESTIMATE SE    p-val    
## nu                 0.281   0.139          
## rho                0.780   0.128          
## lambda             0.010   0.008          
## diabetes_yn_10     1.683   0.970 0.083 .  
## tobacco_use_yn_10  0.117   0.758 0.878    
## male_yn_10        -0.311   0.865 0.719    
## age_stdzd         -0.932   0.554 0.093 .  
## bop_stdzd         -0.074   0.196 0.705    
## plaque_stdzd      -0.382   0.375 0.309    
## calmean_stdzd      0.323   0.241 0.18     
## crown_yn_10        2.024   0.502 <.001 ***
## molar_yn_10        0.049   0.428 0.909    
## ---
## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Kendall's Tau: 0.281
```
## Conditional Proportional Hazards (Stirling-Static in R)


```r
start_time <- Sys.time()
start_time
```

```
## [1] "2020-11-13 15:40:21 EST"
```

```r
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
           integral_term_0007(esa,sa,s,a),
           integral_term_0008(esa,sa,s,a),
           integral_term_0009(esa,sa,s,a),

           integral_term_0010(esa,sa,s,a),
           integral_term_0011(esa,sa,s,a),
           integral_term_0012(esa,sa,s,a),
           integral_term_0013(esa,sa,s,a),
           integral_term_0014(esa,sa,s,a),
           integral_term_0015(esa,sa,s,a),
           integral_term_0016(esa,sa,s,a),
           integral_term_0017(esa,sa,s,a),
           integral_term_0018(esa,sa,s,a),
           integral_term_0019(esa,sa,s,a),

           integral_term_0020(esa,sa,s,a),
           integral_term_0021(esa,sa,s,a),
           integral_term_0022(esa,sa,s,a),
           integral_term_0023(esa,sa,s,a),
           integral_term_0024(esa,sa,s,a),
           integral_term_0025(esa,sa,s,a),
           integral_term_0026(esa,sa,s,a),
           integral_term_0027(esa,sa,s,a),
           integral_term_0028(esa,sa,s,a),
           integral_term_0029(esa,sa,s,a),

           integral_term_0030(esa,sa,s,a),
           integral_term_0031(esa,sa,s,a),
           integral_term_0032(esa,sa,s,a)
           )
}

already_marg_like3 <- function(parms){
  
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
  # beta10   <- parms[13];
  # beta11   <- parms[14];
  # beta12   <- parms[15];

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

## check
options(digits=10)
already_marg_like3(start_parms_val)
```

```
## [1] 205.3008739
```

```r
sprintf("%.10f",sum(round(already_marg_like3(start_parms_val),6)))
```

```
## [1] "205.3008740000"
```

```r
cphz_stirling <-
optim(start_parms_val,
      already_marg_like3,
      method="L-BFGS-B",
      
      lower=c(0.1,0.1, -7, rep(-2, 9)),  ## COND  PH lower bounds     
      upper=c(0.9,2.0,  3, rep( 2, 9)),  ## COND  PH upper bounds     

      hessian=TRUE,
      control=list(factr=10^3.9))
sprintf("%.10f", -cphz_stirling$value)
```

```
## [1] "-100.5687669536"
```

```r
sprintf("%.10f",attr(cphz_recomega, "loglik"))
```

```
## [1] "-100.5676297313"
```

```r
## YES!
fit<-cphz_stirling
fisher_info <- solve(fit$hessian)
## fisher_info <- solve(numDeriv::hessian(already_marg_like3, start_parms_val))
prop_sigma<-matrix(sqrt(diag(fisher_info)))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
already_marg_results <-data.frame(value=fit$par, low95=lower, upp95=upper, loglik=sprintf("%.10f", -cphz_stirling$value))


end_time <- Sys.time()
end_time - start_time
```

```
## Time difference of 3.427079916 secs
```

```r
already_marg_results
```

```
##             value         low95         upp95          loglik
## 1   0.71986189194  0.4467264901  0.9929972938 -100.5687669536
## 2   0.77861422170  0.5276528288  1.0295756146 -100.5687669536
## 3  -4.60660146137 -6.3386601423 -2.8745427804 -100.5687669536
## 4   1.67900411129 -0.2215122959  3.5795205185 -100.5687669536
## 5   0.11438015551 -1.3973475389  1.6261078499 -100.5687669536
## 6  -0.31615173219 -2.0231166199  1.3908131555 -100.5687669536
## 7  -0.92629800188 -2.0113584405  0.1587624367 -100.5687669536
## 8  -0.07229727216 -0.4582825747  0.3136880304 -100.5687669536
## 9  -0.38552047236 -1.1192340618  0.3481931170 -100.5687669536
## 10  0.32181271163 -0.1498420200  0.7934674433 -100.5687669536
## 11  2.00000000000  1.0071528217  2.9928471783 -100.5687669536
## 12  0.05107990299 -0.7905180585  0.8926778645 -100.5687669536
```
