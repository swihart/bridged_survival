proc import out=teeth
	datafile="../diabetic_dental_data.csv"
	dbms=csv replace;
run;

ods html body='conditional_accelerated_failure_200__GQ_gconv0_noadscale_FTOL.htm'(title="cond_aft_200__GQ_gconv0_noadscale_FTOL") style=HTMLBlue;


proc means data=teeth ;
run;

proc nlmixed data=teeth qpoints=200 NTHREADS=-1 noad noadscale gconv=0 FTOL=0.000000000001763762991;

    parms              alpha     =  0.50
                       shape_c   =  1.00
                       gamma_c_f =  5.17
               diabetes_beta_c_f = -0.35
                tobacco_beta_c_f = -1.27
                   male_beta_c_f = -0.20
              age_stdzd_beta_c_f = -0.22
              bop_stdzd_beta_c_f = -0.12
           plaque_stdzd_beta_c_f =  0.02
          calmean_stdzd_beta_c_f = -0.33
                  crown_beta_c_f = -0.18
                  molar_beta_c_f = -0.34;

    bounds               0.10 <              alpha     <   0.90,
                         0.10 <              shape_c   <   2.00,
                         4.92 <              gamma_c_f <   5.42,
                        -0.60 <      diabetes_beta_c_f <  -0.10,
                        -1.52 <       tobacco_beta_c_f <  -1.02,
                        -0.45 <          male_beta_c_f <   0.05,
                        -0.47 <     age_stdzd_beta_c_f <   0.03,
                        -0.37 <     bop_stdzd_beta_c_f <   0.13,
                        -0.23 <  plaque_stdzd_beta_c_f <   0.27,
                        -0.58 < calmean_stdzd_beta_c_f <  -0.08,
                        -0.43 <         crown_beta_c_f <   0.07,
                        -0.59 <         molar_beta_c_f <  -0.09;

    /*setup constants*/
    pi=constant("pi");

    /* generate independent k1~Uniform(0,pi) and k2~Exponential(1)*/
    k1=     probnorm(z1)*pi;
    k2=-log(probnorm(z2));
    /* generate u_i using k1, k2*/
    u_i = sin(alpha * k1) /( sin(k1)**(1/alpha)) * (sin((1-alpha)*k1)/k2)**(1/alpha-1);

    /* linear predictor*/
    lin_pred =                     gamma_c_f                  +
                           diabetes_beta_c_f *    diabetes_yn_10 +
                            tobacco_beta_c_f * tobacco_use_yn_10 +
                               male_beta_c_f *        male_yn_10 +
                          age_stdzd_beta_c_f *         age_stdzd +
                          bop_stdzd_beta_c_f *         bop_stdzd +
                       plaque_stdzd_beta_c_f *      plaque_stdzd +
                      calmean_stdzd_beta_c_f *     calmean_stdzd +
                              crown_beta_c_f *       crown_yn_10 +
                              molar_beta_c_f *       molar_yn_10 ;

    /* parameterize shape_c and scale*/
    scale =  exp(-shape_c * lin_pred);  /* COND AFT parameterization */



    /*likelihood calculations*/
    log_S_t_u_i =  (- u_i * scale * time**shape_c);
    log_h_t_u_i = log(u_i * scale * shape_c * time**(shape_c-1));
    log_f_t_u_i = log_h_t_u_i + log_S_t_u_i;

    if status=0 then loglik = log_S_t_u_i;
    if status=1 then loglik = log_f_t_u_i;

    model time~general(loglik);

    random z1 z2 ~ normal([0,0],[1,0,1]) subject=id;


run;

ods html close;
