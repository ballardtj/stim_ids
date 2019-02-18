functions{

     real lba_pdf(real t, real b, real A, real v, real s){
          //PDF of the LBA model

          real b_A_tv_ts;
          real b_tv_ts;
          real term_1;
          real term_2;
          real term_3;
          real term_4;
          real pdf;

          b_A_tv_ts = (b - A - t*v)/(t*s);
          b_tv_ts = (b - t*v)/(t*s);
          term_1 = v*Phi_approx(b_A_tv_ts);
          term_2 = s*exp(normal_lpdf(b_A_tv_ts|0,1));
          term_3 = v*Phi_approx(b_tv_ts);
          term_4 = s*exp(normal_lpdf(b_tv_ts|0,1));
          pdf = (1/A)*(-term_1 + term_2 + term_3 - term_4);

          return pdf;
     }

     real lba_cdf(real t, real b, real A, real v, real s){
          //CDF of the LBA model

          real b_A_tv;
          real b_tv;
          real ts;
          real term_1;
          real term_2;
          real term_3;
          real term_4;
          real cdf;

          b_A_tv = b - A - t*v;
          b_tv = b - t*v;
          ts = t*s;
          term_1 = b_A_tv/A * Phi_approx(b_A_tv/ts);
          term_2 = b_tv/A   * Phi_approx(b_tv/ts);
          term_3 = ts/A     * exp(normal_lpdf(b_A_tv/ts|0,1));
          term_4 = ts/A     * exp(normal_lpdf(b_tv/ts|0,1));
          cdf = 1 + term_1 - term_2 + term_3 - term_4;

          return cdf;

     }

     real lba_log(matrix RT, //LENGTH X 2 matrix of RT and response data
                  vector A,    //start point variability parameter for each subject
                  vector B,    //threshold parameter for each subject
                  real v_false, //mean drift for false accumulator for each subject
                  vector s,    //sd of drift rate parameter
                  vector tau,
                  vector v_true){

          real t;
          vector[rows(RT)] b;
          //vector b;
          real cdf;
          real pdf;
          vector[rows(RT)] prob;
          real out;
          real prob_neg;
          real v[2];

          v[1] = v_false;
          b = A + B;

          for (i in 1:rows(RT)){
              v[2] = v_true[i];
               t = RT[i,1] - tau[i];
               if(t > 0){
                    cdf = 1;
                    for(j in 1:num_elements(v)){
                         if(RT[i,2] == j){
                              pdf = lba_pdf(t, b[i], A[i], v[j], s[i]);
                         }else{
                              cdf = (1-lba_cdf(t, b[i], A[i], v[j], s[i])) * cdf;
                         }
                    }
                    prob_neg = 1;
                    for(j in 1:num_elements(v)){
                         prob_neg = Phi_approx(-v[j]/s[i]) * prob_neg;
                    }
                    prob[i] = pdf*cdf;
                    prob[i] = prob[i]/(1-prob_neg);
                    if(prob[i] < 1e-10){
                         prob[i] = 1e-10;
                    }

               }else{
                    prob[i] = 1e-10;
               }
          }
          out = sum(log(prob));
          return out;
     }

    vector lba_rng(real k, real A, vector v, real s, real tau){

          int get_pos_drift;
          int no_pos_drift;
          int get_first_pos;
          vector[num_elements(v)] drift;
          int max_iter;
          int iter;
          real start[num_elements(v)];
          real ttf[num_elements(v)];
          int resp[num_elements(v)];
          real rt;
          vector[2] pred;
          real b;

          //try to get a positive drift rate
          get_pos_drift = 1;
          no_pos_drift = 0;
          max_iter = 1000;
          iter = 0;
          while(get_pos_drift){
               for(j in 1:num_elements(v)){
                    drift[j] = normal_rng(v[j],s);
                    if(drift[j] > 0){
                         get_pos_drift = 0;
                    }
               }
               iter = iter + 1;
               if(iter > max_iter){
                    get_pos_drift = 0;
                    no_pos_drift = 1;
               }
          }
          //if both drift rates are <= 0
          //return an infinite response time
          if(no_pos_drift){
               pred[1] = -1;
               pred[2] = -1;
          }else{
               b = A + k;
               for(i in 1:num_elements(v)){
                    //start time of each accumulator
                    start[i] = uniform_rng(0,A);
                    //finish times
                    ttf[i] = (b-start[i])/drift[i];
               }
               //rt is the fastest accumulator finish time
               //if one is negative get the positive drift
               resp = sort_indices_asc(ttf);
               ttf = sort_asc(ttf);
               get_first_pos = 1;
               iter = 1;
               while(get_first_pos){
                    if(ttf[iter] > 0){
                         pred[1] = ttf[iter] + tau;
                         pred[2] = resp[iter];
                         get_first_pos = 0;
                    }
                    iter = iter + 1;
               }
          }
          return pred;
     }
}

data{
     int Ntotal;
     int Nsubj;
     int Nvars;
     int subject[Ntotal];
     matrix[Nsubj,Nvars] COVS;
     vector[Ntotal] dpost_anodal;
     vector[Ntotal] dpost_cathodal;
     vector[Ntotal] dpost_sham;
     //
     // matrix[Ntotal,Nvars] dpost_anodal_vars;
     // matrix[Ntotal,Nvars] dpost_cathodal_vars;
     // matrix[Ntotal,Nvars] dpost_sham_vars;
     matrix[Ntotal,2] RT;
}

// transformed data{
//   vector[Ntotal] dpost_anodal_on_2;
//   vector[Ntotal] dpost_cathodal_on_2;
//   vector[Ntotal] dpost_sham_on_2;
//
//   dpost_anodal_on_2 = dpost_anodal / 2;
//   dpost_cathodal_on_2 = dpost_cathodal / 2;
//   dpost_sham_on_2 = dpost_sham / 2;
//
// }

parameters {
    //hyperparameters
    real<lower=0> A_mean;
    real<lower=0> A_sd;
    real<lower=0> s_mean;
    real<lower=0> s_sd;
    real<lower=0> B_pre_mean;
    real<lower=0> B_pre_sd;
    real<lower=0.1,upper=1> tau_mean;
    real<lower=0,upper=1> tau_sd;

    real v_true_pre_mean;
    real<lower=0> v_true_pre_sd;

    real dvt_sham_mean;
    real<lower=0> dvt_sham_sd;
    real dB_sham_mean;
    real<lower=0> dB_sham_sd;

    vector[Nvars] COEFS_anodal_vt;
    real<lower=0> vt_anodal_sd;
    vector[Nvars] COEFS_cathodal_vt;
    real<lower=0> vt_cathodal_sd;

    vector[Nvars] COEFS_anodal_B;
    real<lower=0> B_anodal_sd;
    vector[Nvars] COEFS_cathodal_B;
    real<lower=0> B_cathodal_sd;

    //person-level parameters
    vector<lower=0>[Nsubj] A_raw;
    vector<lower=0>[Nsubj] s_raw;
    vector<lower=0>[Nsubj] B_pre_raw;
    vector<lower=0.1,upper=1>[Nsubj] tau;
    vector[Nsubj] v_true_pre_raw;
    vector[Nsubj] vt_anodal_disruption_raw;   //change in anodal
    vector[Nsubj] vt_cathodal_disruption_raw; //change in cathodal condition
    vector[Nsubj] dvt_sham_raw;     //change in sham condition
    vector[Nsubj] B_anodal_disruption_raw;   //change in anodal
    vector[Nsubj] B_cathodal_disruption_raw; //change in cathodal condition
    vector[Nsubj] dB_sham_raw;     //change in sham condition

}

transformed parameters {
     vector[Nsubj] s;
     vector<lower=0>[Nsubj] A;
     //vector[Nsubj] v_false_pre;
     vector[Nsubj] v_true_pre;
     vector<lower=0>[Nsubj] B_pre;

     //Outputs of regression models
     vector[Nsubj] vt_anodal_disruption_mean;
     vector[Nsubj] vt_cathodal_disruption_mean;
     vector[Nsubj] B_anodal_disruption_mean;
     vector[Nsubj] B_cathodal_disruption_mean;

     vector[Nsubj] vt_anodal_disruption;
     vector[Nsubj] vt_cathodal_disruption;
     vector[Nsubj] B_anodal_disruption;
     vector[Nsubj] B_cathodal_disruption;

     vector[Nsubj] dvt_anodal;
     vector[Nsubj] dvt_cathodal;
     vector[Nsubj] dvt_sham;

     vector[Nsubj] dB_anodal;
     vector[Nsubj] dB_cathodal;
     vector[Nsubj] dB_sham;

     vector[Ntotal] dB;
     real dB_min;

     vector[Ntotal] v_true;
     vector[Ntotal] B;

     //set constant
     real v_false = 1;

     //uncenter parameters
     A = A_mean + A_sd*A_raw;
     s = s_mean + s_sd*s_raw;
     B_pre = B_pre_mean + B_pre_sd*B_pre_raw;
     v_true_pre = v_true_pre_mean + v_true_pre_sd*v_true_pre_raw;
     dvt_sham = dvt_sham_mean + dvt_sham_sd*dvt_sham_raw;
     dB_sham = dB_sham_mean + dB_sham_sd*dB_sham_raw;


     //calculate means for distribution parameters
     vt_anodal_disruption_mean = COVS * COEFS_anodal_vt;
     vt_cathodal_disruption_mean = COVS * COEFS_cathodal_vt;
     B_anodal_disruption_mean = COVS * COEFS_anodal_B;
     B_cathodal_disruption_mean = COVS * COEFS_cathodal_B;

     //uncenter distruption variables
     vt_anodal_disruption = vt_anodal_disruption_raw*vt_anodal_sd + vt_anodal_disruption_mean;
     vt_cathodal_disruption = vt_cathodal_disruption_raw*vt_cathodal_sd + vt_cathodal_disruption_mean;
     B_anodal_disruption = B_anodal_disruption_raw*B_anodal_sd + B_anodal_disruption_mean;
     B_cathodal_disruption = B_cathodal_disruption_raw*B_cathodal_sd + B_cathodal_disruption_mean;

     //calculate vt and B post for anodal and cathodal.
     dvt_anodal = dvt_sham + vt_anodal_disruption;
     dvt_cathodal = dvt_sham + vt_cathodal_disruption;
     dB_anodal = dB_sham + B_anodal_disruption;
     dB_cathodal = dB_sham + B_cathodal_disruption;

    //calculate v_true
    v_true = v_true_pre[subject] +
        dvt_anodal[subject] .* dpost_anodal +
        dvt_cathodal[subject] .* dpost_cathodal +
        dvt_sham[subject] .* dpost_sham;

    //calculate B
    dB = dB_anodal[subject] .* dpost_anodal +
        dB_cathodal[subject] .* dpost_cathodal +
        dB_sham[subject] .* dpost_sham;

     dB_min = min(dB); //calculates lowest threshold

     B = B_pre[subject] - dB_min + dB; //actual B_pre = B_pre - min; function is "pushed up" so threshold must be positive

}

model {
     //hyperpriors
     A_mean ~ normal(1,1); //A_sd ~ uniform(0,1)
     A_sd ~ normal(0,1);
     s_mean ~ normal(1,1);
     s_sd ~ normal(0,1);

     B_pre_mean ~ normal(1,1); //B_pre_sd ~ uniform(0,1)
     B_pre_sd ~ normal(0,1);
     v_true_pre_mean ~ normal(1,2); //v_true_pre_sd ~ uniform(0,1)
     v_true_pre_sd ~ normal(0,1);

     dvt_sham_mean ~ normal(0,1);
     dvt_sham_sd ~ normal(0,1);
     dB_sham_mean ~ normal(0,1);
     dB_sham_sd ~ normal(0,1);


     COEFS_anodal_vt ~ normal(0,2);
     COEFS_cathodal_vt ~ normal(0,2);
     vt_anodal_sd ~ normal(0,1);
     vt_cathodal_sd ~ normal(0,1);

     COEFS_anodal_B ~ normal(0,2);
     COEFS_cathodal_B ~ normal(0,2);
     B_anodal_sd ~ normal(0,1);
     B_cathodal_sd ~ normal(0,1);

     //priors
     A_raw ~ normal(0,1);
     s_raw ~ normal(0,1);
     B_pre_raw ~ normal(0,1);
    // v_false_pre_raw ~ normal(0,1);
     v_true_pre_raw ~ normal(0,1);
     //tau_sd ~ normal(0,1);

     for(subj in 1:Nsubj){
       tau[subj] ~ normal(tau_mean,tau_sd)T[0.1,1];
     }

     //regression (not sure if these can be uncentred)
     vt_anodal_disruption_raw ~ normal(0,1);  //  implies  normal(dvt_anodal_mean,dvt_anodal_sd);
     vt_cathodal_disruption_raw ~ normal(0,1);  //  implies normal(dvt_cathodal_mean,dvt_cathodal_sd);
     dvt_sham_raw ~ normal(0,1);  //  implies normal(dvt_sham_mean,dvt_sham_sd);

     B_anodal_disruption_raw ~ normal(0,1);  //  implies normal(dB_anodal_mean,dB_anodal_sd);
     B_cathodal_disruption_raw ~ normal(0,1);  //  implies normal(dB_cathodal_mean,dB_cathodal_sd);
     dB_sham_raw ~ normal(0,1);  //  implies normal(dB_sham_mean,dB_sham_sd);

     RT ~ lba(A[subject],B,v_false,s[subject],tau[subject],v_true);
}






