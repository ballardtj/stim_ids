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
                  vector v_false, //mean drift for false accumulator for each subject
                  real s,    //sd of drift rate parameter
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
          //v[1] = v_false;
          //v[2] = v_true;
          b = A + B;

          for (i in 1:rows(RT)){
              v[1] = v_false[i];
              v[2] = v_true[i];
               t = RT[i,1] - tau[i];
               if(t > 0){
                    cdf = 1;
                    for(j in 1:num_elements(v)){
                         if(RT[i,2] == j){
                              pdf = lba_pdf(t, b[i], A[i], v[j], s);
                         }else{
                              cdf = (1-lba_cdf(t, b[i], A[i], v[j], s)) * cdf;
                         }
                    }
                    prob_neg = 1;
                    for(j in 1:num_elements(v)){
                         prob_neg = Phi_approx(-v[j]/s) * prob_neg;
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
    //real<lower=0,upper=1> A_sd;
    real<lower=0> A_sd;
    real<lower=0> B_pre_mean;
    real<lower=0> B_pre_sd;
    //real<lower=0,upper=1> B_pre_sd;
    real<lower=0.1,upper=1> tau_mean;
    real<lower=0> tau_sd;
    //real<lower=0,upper=1> tau_sd;
    real v_false_pre_mean;
    real<lower=0> v_false_pre_sd;
    //real<lower=0,upper=1> v_false_pre_sd;
    real v_true_pre_mean;
    real<lower=0> v_true_pre_sd;
    //real<lower=0,upper=1> v_true_pre_sd;

    vector[Nvars] COEFS_anodal_v_true;
    real<lower=0> dvt_anodal_sd;
    vector[Nvars] COEFS_cathodal_v_true;
    real<lower=0> dvt_cathodal_sd;
    vector[Nvars] COEFS_sham_v_true;
    real<lower=0> dvt_sham_sd;

    vector[Nvars] COEFS_anodal_v_false;
    real<lower=0> dvf_anodal_sd;
    vector[Nvars] COEFS_cathodal_v_false;
    real<lower=0> dvf_cathodal_sd;
    vector[Nvars] COEFS_sham_v_false;
    real<lower=0> dvf_sham_sd;

    vector[Nvars] COEFS_anodal_B;
    real<lower=0> dB_anodal_sd;
    vector[Nvars] COEFS_cathodal_B;
    real<lower=0> dB_cathodal_sd;
    vector[Nvars] COEFS_sham_B;
    real<lower=0> dB_sham_sd;

    //person-level parameters
    vector<lower=0>[Nsubj] A_raw;
    vector<lower=0>[Nsubj] B_pre_raw;
    vector<lower=0.1,upper=1>[Nsubj] tau;
    vector[Nsubj] v_false_pre_raw;
    vector[Nsubj] v_true_pre_raw;

    vector[Nsubj] dvt_anodal_raw;   //change in anodal
    vector[Nsubj] dvt_cathodal_raw; //change in cathodal condition
    vector[Nsubj] dvt_sham_raw;     //change in sham condition
    vector[Nsubj] dvf_anodal_raw;   //change in anodal
    vector[Nsubj] dvf_cathodal_raw; //change in cathodal condition
    vector[Nsubj] dvf_sham_raw;     //change in sham condition
    vector[Nsubj] dB_anodal_raw;   //change in anodal
    vector[Nsubj] dB_cathodal_raw; //change in cathodal condition
    vector[Nsubj] dB_sham_raw;     //change in sham condition


}

transformed parameters {
     real s;
     vector<lower=0>[Nsubj] A;
     vector[Nsubj] v_false_pre;
     vector[Nsubj] v_true_pre;
     vector<lower=0>[Nsubj] B_pre;
     vector[Nsubj] dvt_anodal_mean;
     vector[Nsubj] dvt_cathodal_mean;
     vector[Nsubj] dvt_sham_mean;
     vector[Nsubj] dvf_anodal_mean;
     vector[Nsubj] dvf_cathodal_mean;
     vector[Nsubj] dvf_sham_mean;
     vector[Nsubj] dB_anodal_mean;
     vector[Nsubj] dB_cathodal_mean;
     vector[Nsubj] dB_sham_mean;
     vector[Nsubj] dvt_anodal;
     vector[Nsubj] dvt_cathodal;
     vector[Nsubj] dvt_sham;
     vector[Nsubj] dvf_anodal;
     vector[Nsubj] dvf_cathodal;
     vector[Nsubj] dvf_sham;
     vector[Nsubj] dB_anodal;
     vector[Nsubj] dB_cathodal;
     vector[Nsubj] dB_sham;
     vector[Ntotal] dB;
     real dB_min;
     vector[Ntotal] v_false;
     vector[Ntotal] v_true;
     vector[Ntotal] B;

     //set constant
     s = 1;

     //uncenter parameters
     A = A_mean + A_sd*A_raw;
     B_pre = B_pre_mean + B_pre_sd*B_pre_raw;
     v_false_pre = v_false_pre_mean + v_false_pre_sd*v_false_pre_raw;
     v_true_pre = v_true_pre_mean + v_true_pre_sd*v_true_pre_raw;

     //calculate means for change parameters
     dvt_anodal_mean = COVS * COEFS_anodal_v_true;
     dvt_cathodal_mean = COVS * COEFS_cathodal_v_true;
     dvt_sham_mean = COVS * COEFS_sham_v_true;

     dvf_anodal_mean = COVS * COEFS_anodal_v_false;
     dvf_cathodal_mean = COVS * COEFS_cathodal_v_false;
     dvf_sham_mean = COVS * COEFS_sham_v_false;

     dB_anodal_mean = COVS * COEFS_anodal_B;
     dB_cathodal_mean = COVS * COEFS_cathodal_B;
     dB_sham_mean = COVS * COEFS_sham_B;

     dvt_anodal = dvt_anodal_raw*dvt_anodal_sd + dvt_anodal_mean;
     dvt_cathodal = dvt_cathodal_raw*dvt_cathodal_sd + dvt_cathodal_mean;
     dvt_sham = dvt_sham_raw*dvt_sham_sd + dvt_sham_mean;

     dvf_anodal = dvf_anodal_raw*dvf_anodal_sd + dvf_anodal_mean;
     dvf_cathodal = dvf_cathodal_raw*dvf_cathodal_sd + dvf_cathodal_mean;
     dvf_sham = dvf_sham_raw*dvf_sham_sd + dvf_sham_mean;

     dB_anodal = dB_anodal_raw*dB_anodal_sd + dB_anodal_mean;
     dB_cathodal = dB_cathodal_raw*dB_cathodal_sd + dB_cathodal_mean;
     dB_sham = dB_sham_raw*dB_sham_sd + dB_sham_mean;



     //calculate v_true
     v_true = v_true_pre[subject] +
        dvt_anodal[subject] .* dpost_anodal +
        dvt_cathodal[subject] .* dpost_cathodal +
        dvt_sham[subject] .* dpost_sham;

    //calculate v_false
     v_false = v_false_pre[subject] +
        dvf_anodal[subject] .* dpost_anodal +
        dvf_cathodal[subject] .* dpost_cathodal +
        dvf_sham[subject] .* dpost_sham;

    //calculate B
    dB = dB_anodal[subject] .* dpost_anodal +
        dB_cathodal[subject] .* dpost_cathodal +
        dB_sham[subject] .* dpost_sham;

     dB_min = min(dB); //calculates lowest change value

     B = B_pre[subject] - dB_min + dB; //actual B_pre = B_pre - min; function is "pushed up" so threshold must be positive




//
//      dvt_anodal = dpost_anodal_vars * COEFS_anodal;
//      dvt_cathodal = dpost_cathodal_vars * COEFS_cathodal;
//      dvt_sham = dpost_sham_vars * COEFS_sham;
//
//      v_true = v_true_pre[subject] + dvt_anodal + dvt_cathodal + dvt_sham;
}

model {
     //hyperpriors
     A_mean ~ normal(1,1); //A_sd ~ uniform(0,1)
     A_sd ~ normal(0,1);
     B_pre_mean ~ normal(1,1); //B_pre_sd ~ uniform(0,1)
     B_pre_sd ~ normal(0,1);
     v_true_pre_mean ~ normal(1,2); //v_true_pre_sd ~ uniform(0,1)
     v_true_pre_sd ~ normal(0,1);
     v_false_pre_mean ~ normal(1,2); //v_false_pre_sd ~ uniform(0,1)
     v_false_pre_sd ~ normal(0,1);
     //v_false_mean ~ normal(1,2); //v_false_sd ~ uniform(0,1)
     //tau_mean ~ uniform(0.1,1); tau_sd ~ uniform(0,1)
     COEFS_anodal_v_true ~ normal(0,2);
     COEFS_cathodal_v_true ~ normal(0,2);
     COEFS_sham_v_true ~ normal(0,2);
     COEFS_anodal_v_false ~ normal(0,2);
     COEFS_cathodal_v_false ~ normal(0,2);
     COEFS_sham_v_false ~ normal(0,2);
     COEFS_anodal_B ~ normal(0,2);
     COEFS_cathodal_B ~ normal(0,2);
     COEFS_sham_B ~ normal(0,2);
     dvt_anodal_sd ~ normal(0,1);
     dvt_cathodal_sd ~ normal(0,1);
     dvt_sham_sd ~ normal(0,1);
     dvf_anodal_sd ~ normal(0,1);
     dvf_cathodal_sd ~ normal(0,1);
     dvf_sham_sd ~ normal(0,1);
     dB_anodal_sd ~ normal(0,1);
     dB_cathodal_sd ~ normal(0,1);
     dB_sham_sd ~ normal(0,1);

     //priors
     A_raw ~ normal(0,1);
     B_pre_raw ~ normal(0,1);
     v_false_pre_raw ~ normal(0,1);
     v_true_pre_raw ~ normal(0,1);
     tau_sd ~ normal(0,1);

     for(subj in 1:Nsubj){
       tau[subj] ~ normal(tau_mean,tau_sd)T[0.1,1];
     }

     //regression (not sure if these can be uncentred)
     dvt_anodal_raw ~ normal(0,1);  //  implies  normal(dvt_anodal_mean,dvt_anodal_sd);
     dvt_cathodal_raw ~ normal(0,1);  //  implies normal(dvt_cathodal_mean,dvt_cathodal_sd);
     dvt_sham_raw ~ normal(0,1);  //  implies normal(dvt_sham_mean,dvt_sham_sd);

     dvf_anodal_raw ~ normal(0,1);  //  implies normal(dvf_anodal_mean,dvf_anodal_sd);
     dvf_cathodal_raw ~ normal(0,1);  //  implies normal(dvf_cathodal_mean,dvf_cathodal_sd);
     dvf_sham_raw ~ normal(0,1);  //  implies normal(dvf_sham_mean,dvf_sham_sd);

     dB_anodal_raw ~ normal(0,1);  //  implies normal(dB_anodal_mean,dB_anodal_sd);
     dB_cathodal_raw ~ normal(0,1);  //  implies normal(dB_cathodal_mean,dB_cathodal_sd);
     dB_sham_raw ~ normal(0,1);  //  implies normal(dB_sham_mean,dB_sham_sd);

     RT ~ lba(A[subject],B,v_false,s,tau[subject],v_true);
}






