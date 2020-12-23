proc import out=teeth
	datafile="../diabetic_dental_data.csv"
	dbms=csv replace;
run;

ods html body='conditional_proportional_hazards_200__GQ_gconv0_noadscale_FTOL.htm'(title="cond__ph_200__GQ_gconv0_noadscale_FTOL") style=HTMLBlue;


proc means data=teeth ;
run;

proc nlmixed data=teeth qpoints=200 NTHREADS=-1 noad noadscale gconv=0 FTOL=0.000000000001763762991;

    parms              alpha     =  0.50
                       shape_c   =  1.00
                       gamma_c_z =  0.10
               diabetes_beta_c_z =  0.10
                tobacco_beta_c_z =  0.10
                   male_beta_c_z =  0.10
              age_stdzd_beta_c_z =  0.10
              bop_stdzd_beta_c_z =  0.10
           plaque_stdzd_beta_c_z =  0.10
          calmean_stdzd_beta_c_z =  0.10
                  crown_beta_c_z =  0.10
                  molar_beta_c_z =  0.10;

    bounds             0.10 <              alpha     < 0.90,
                       0.10 <              shape_c   < 2.00,
                      -6.00 <              gamma_c_z < 6.00,
                      -2.00 <      diabetes_beta_c_z < 2.00,
                      -2.00 <       tobacco_beta_c_z < 2.00,
                      -2.00 <          male_beta_c_z < 2.00,
                      -2.00 <     age_stdzd_beta_c_z < 2.00,
                      -2.00 <     bop_stdzd_beta_c_z < 2.00,
                      -2.00 <  plaque_stdzd_beta_c_z < 2.00,
                      -2.00 < calmean_stdzd_beta_c_z < 2.00,
                      -2.00 <         crown_beta_c_z < 2.00,
                      -2.00 <         molar_beta_c_z < 2.00;


    /*setup constants*/
    pi=constant("pi");

    /* generate independent k1~Uniform(0,pi) and k2~Exponential(1)*/
    k1=     probnorm(z1)*pi;
    k2=-log(probnorm(z2));
    /* generate u_i using k1, k2*/
    u_i = sin(alpha * k1) /( sin(k1)**(1/alpha)) * (sin((1-alpha)*k1)/k2)**(1/alpha-1);

    /* linear predictor*/
    lin_pred =                     gamma_c_z                  +
                           diabetes_beta_c_z *    diabetes_yn_10 +
                            tobacco_beta_c_z * tobacco_use_yn_10 +
                               male_beta_c_z *        male_yn_10 +
                          age_stdzd_beta_c_z *         age_stdzd +
                          bop_stdzd_beta_c_z *         bop_stdzd +
                       plaque_stdzd_beta_c_z *      plaque_stdzd +
                      calmean_stdzd_beta_c_z *     calmean_stdzd +
                              crown_beta_c_z *       crown_yn_10 +
                              molar_beta_c_z *       molar_yn_10 ;

    /* parameterize shape_c and scale*/
    scale =  exp(lin_pred);         /* COND PH  parameterization */


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
