proc import out=teeth
	datafile="../diabetic_dental_data.csv"
	dbms=csv replace;
run;

ods html body='marginal_proportional_hazards_200__GQ_gconv0_noadscale_FTOL.htm'(title="marg__ph_200__GQ_gconv0_noadscale_FTOL") style=HTMLBlue;


proc means data=teeth ;
run;

proc nlmixed data=teeth qpoints=200 NTHREADS=-1 noad noadscale gconv=0 FTOL=0.000000000001763762991;

    parms              alpha     =  0.50
                       shape_m   =  1.00
                       gamma_m_z = -3.22
               diabetes_beta_m_z =  0.22
                tobacco_beta_m_z =  0.79
                   male_beta_m_z =  0.12
              age_stdzd_beta_m_z =  0.14
              bop_stdzd_beta_m_z =  0.07
           plaque_stdzd_beta_m_z = -0.01
          calmean_stdzd_beta_m_z =  0.21
                  crown_beta_m_z =  0.11
                  molar_beta_m_z =  0.21;

    bounds             0.1 <              alpha     <  0.9,
                       0.1 <              shape_m   <  2.0,
                      -4.0 <              gamma_m_z < -3.0,
                       0.0 <      diabetes_beta_m_z <  1.0,
                       0.0 <       tobacco_beta_m_z <  1.0,
                       0.0 <          male_beta_m_z <  1.0,
                       0.0 <     age_stdzd_beta_m_z <  1.0,
                       0.0 <     bop_stdzd_beta_m_z <  1.0,
                      -1.0 <  plaque_stdzd_beta_m_z <  0.0,
                       0.0 < calmean_stdzd_beta_m_z <  1.0,
                       0.0 <         crown_beta_m_z <  1.0,
                       0.0 <         molar_beta_m_z <  1.0;

    /*setup constants*/
    pi=constant("pi");

    /* generate independent k1~Uniform(0,pi) and k2~Exponential(1)*/
    k1=     probnorm(z1)*pi;
    k2=-log(probnorm(z2));
    /* generate u_i using k1, k2*/
    u_i = sin(alpha * k1) /( sin(k1)**(1/alpha)) * (sin((1-alpha)*k1)/k2)**(1/alpha-1);

    /* linear predictor*/
    lin_pred =                     gamma_m_z                  +
                           diabetes_beta_m_z *    diabetes_yn_10 +
                            tobacco_beta_m_z * tobacco_use_yn_10 +
                               male_beta_m_z *        male_yn_10 +
                          age_stdzd_beta_m_z *         age_stdzd +
                          bop_stdzd_beta_m_z *         bop_stdzd +
                       plaque_stdzd_beta_m_z *      plaque_stdzd +
                      calmean_stdzd_beta_m_z *     calmean_stdzd +
                              crown_beta_m_z *       crown_yn_10 +
                              molar_beta_m_z *       molar_yn_10 ;

    /* parameterize shape_c and scale*/
    shape_c = shape_m / alpha;
    scale =  exp(lin_pred/alpha);  /* MARG PH parameterization */



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
