data {
  //number of sites
  int<lower = 1> n_site;
  int<lower = 1> n_year;
  int<lower = 1> n_mtb;
  int<lower = 1> n_state;
  
  // survey-level detection covariates - multiple surveys at each site
  int<lower = 1> total_surveys;
  int<lower = 1> m_p;// number of covariates
  matrix[total_surveys, m_p] X_p;// detection covariate matrix
  int detYear[total_surveys];
  int detState[total_surveys];
  
  // survey level information 
  int<lower = 1, upper = n_site> site[total_surveys];
  int<lower = 0, upper = 1> y[total_surveys];//survey-level observation - presence/absence
  int<lower = 0, upper = total_surveys> start_idx[n_site];
  int<lower = 0, upper = total_surveys> end_idx[n_site];
  int<lower = 0, upper = total_surveys> start_midx[n_mtb];
  int<lower = 0, upper = total_surveys> end_midx[n_mtb];

  // spline data for occupancy model for site variation
  // data for covariates
  int m_psi;  // number of covariates
  matrix[n_site, m_psi] X_psi;  // design matrix for the covariates
  // data for splines
  int Ks;  // number of linear effects
  matrix[n_site, Ks] Xs;  // design matrix for the linear effects
  // data for spline t2(x,y)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[n_site, knots_1[1]] Zs_1_1;
  matrix[n_site, knots_1[2]] Zs_1_2;
  matrix[n_site, knots_1[3]] Zs_1_3;
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[n_site];
  int<lower = 0> total_seen[n_site];
  int<lower = 0> total_seen_mtb[n_mtb];
  
  // number of surveys at each site
  int<lower = 0> n_survey[n_site];
  
  //data for predictions
  int<lower = 1> complete_N;
  // data for covariates
  matrix[complete_N, m_psi] X_psi_complete; 
  // data for splines
  int complete_Ks;  // number of linear effects
  matrix[complete_N, complete_Ks] complete_Xs;  // design matrix for the linear effects
  // data for spline t2(x,y,yearIndex)
  int complete_nb_1;  // number of bases
  int complete_knots_1[complete_nb_1];  // number of knots
  // basis function matrices
  matrix[complete_N, complete_knots_1[1]] complete_Zs_1_1;
  matrix[complete_N, complete_knots_1[2]] complete_Zs_1_2;
  matrix[complete_N, complete_knots_1[3]] complete_Zs_1_3;
}  
parameters {
  
  //detection covariates
  vector[m_p] beta_p;
  real year[n_year];
  real state[n_state];
  real<lower=0,upper=10> yearsigma;
  real<lower=0,upper=10> statesigma;
  
  //occcupancy covariates
  vector[m_psi] beta_psi;
  
  //spline parameters
  vector[Ks] bs;  // spline coefficients
  // parameters for spline t2(x,y,yearIndex)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  vector[knots_1[2]] zs_1_2;
  vector[knots_1[3]] zs_1_3;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  real<lower=0> sds_1_2;  // standard deviations of spline coefficients
  real<lower=0> sds_1_3;  // standard deviations of spline coefficients
  
}
transformed parameters {
  
  real<lower=0,upper=1>  detp[total_surveys];
  real logit_p[total_surveys];
  
  //spline transformed parameters
  vector[knots_1[1]] s_1_1;
  vector[knots_1[2]] s_1_2;
  vector[knots_1[3]] s_1_3;
  s_1_1 = sds_1_1 * zs_1_1;
  s_1_2 = sds_1_2 * zs_1_2;
  s_1_3 = sds_1_3 * zs_1_3;
  
  // detection linear predictor
  for(i in 1:total_surveys){
    detp[i] = inv_logit(X_p[i,] * beta_p + year[detYear[i]] + state[detState[i]]);
  };
  logit_p = logit(detp);
  
}
model {
  
  //occupancy model linear predictor
  vector[n_site] logit_psi = X_psi * beta_psi + Xs * bs + Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3;
  
  //transform to log scale
  vector[n_site] log_psi = log_inv_logit(logit_psi);
  vector[n_site] log1m_psi = log1m_inv_logit(logit_psi);
  
  //priors for detection model
  beta_p ~ normal(0, 2);
  year ~ normal(0, yearsigma);
  state ~ normal(0, statesigma);
  
  //priors for occupancy model
  beta_psi ~ normal(0, 2);
  
  //priors for spline occupancy model
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_1_2);
  target += std_normal_lpdf(zs_1_3);
  
  //likelihood
  for (i in 1:n_site) {
      if (any_seen[i]==1) {
        // if species is seen at the site at least once - site is occupied
        target += log_psi[i] 
        + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] | 
                                 logit_p[start_idx[i]:end_idx[i]]);
      } else {
        // if species never seen at the site - site may or may not be occupied
        target += log_sum_exp(
          log_psi[i] + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] |
                                              logit_p[start_idx[i]:end_idx[i]]), 
          log1m_psi[i]
        );
      }
    }
}
// to change eventually
generated quantities {
  vector[n_site] psi;
  vector[n_site] z;
  vector[n_site] z_rep;
  vector[n_site] top;
  vector[n_site] bottom;
  vector[n_mtb] preds_rep;
  vector[n_mtb] preds_exp;
  vector[n_mtb] sim_rep;
  vector[n_mtb] sim_exp;
  vector[complete_N] psi_comp;
  vector[complete_N] psi_comp_rep;
  vector[total_surveys] p;
  vector[total_surveys] y_rep;
  vector[total_surveys] y_exp;
  real bpv;
  real tot_exp;
  real tot_rep;
  real mean_p;
  real prop_zero;
  
  //psi on original scale for all sites
  for(i in 1:complete_N) { 
    psi_comp[i] = inv_logit(X_psi_complete[i] * beta_psi + complete_Xs[i] * bs + complete_Zs_1_1[i] * s_1_1 + complete_Zs_1_2[i] * s_1_2 + complete_Zs_1_3[i] * s_1_3);
    psi_comp_rep[i] = bernoulli_rng(psi_comp[i]);
  }
  
  // p on original scale
  for(i in 1:total_surveys){
    p[i] = inv_logit(X_p[i,] * beta_p + year[detYear[i]] + state[detState[i]]);
  }
  mean_p = mean(p);
  
  // calc z
  for(i in 1:n_site){
    psi[i] = inv_logit(X_psi[i] * beta_psi + Xs[i] * bs + Zs_1_1[i] * s_1_1 + Zs_1_2[i] * s_1_2 + Zs_1_3[i] * s_1_3);
  
  if(any_seen[i]==1) {
      z[i] = 1; 
  }else {
      top[i] = psi[i] * prod(1- p[start_idx[i]:end_idx[i]]); 
      bottom[i] = (1-psi[i]) + psi[i] * prod(1- p[start_idx[i]:end_idx[i]]);
      z[i] = top[i]/bottom[i];
    }
    z_rep[i] = bernoulli_rng(psi[i]);
  }
  
  for(i in 1:total_surveys){
  y_exp[i] =  z[site[i]] * p[i];
  y_rep[i] =  bernoulli_rng(z_rep[site[i]] * p[i]);
}

//sum per mtb
for(i in 1:n_mtb) {
    preds_rep[i] = sum(y_rep[start_midx[i]:end_midx[i]]);
    preds_exp[i] = sum(y_exp[start_midx[i]:end_midx[i]]);
  }
  
// freeman tukey
for(i in 1:n_mtb){
  sim_exp[i] = square(total_seen_mtb[i] - preds_exp[i]+0.001)/(preds_exp[i]+0.001);
  sim_rep[i] = square(preds_rep[i] - preds_exp[i]+0.001)/(preds_exp[i]+0.001);
}

}
