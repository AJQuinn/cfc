function [results] = estimate_cfc_metrics(signal)


results = [];

%% Estimate the PLV - using theta phase and gamma amplitude
% Lachaux 1999
results.plv = cfc.est.plv(signal.theta_phase, signal.gamma_amp_phase);

%% Estimate AEC - using theta amplitude and gamma amplitude 
% Bruns & Eckhorn 2004
results.aec= cfc_est_aec(signal.theta_amp, signal.gamma_amp);

%% Estimate MI - using theta phase and gamma amplitude
% Canolty 2006
results.mi = cfc.est.mi(signal.theta_phase,signal.gamma_amp);

%% Estimate ESC - using theta signal and gamma amplitude
% Bruns & Eckhorn 2004
results.esc = cfc_est_esc(signal.theta,signal.gamma_amp);

%% Estimate NESC - using theta phase and gamma amplitude
% Penny 2007
results.nesc = ncfc_est_esc(signal.theta_phase, signal.gamma_amp);

%% Estimate GLM -  using theta phase and gamma amplitude
% Penny 2007
results.glm = cfc.est.glm(signal.theta_phase, signal.gamma_amp);

%% Estimate Voytek - using theta phase and theta band filtered gamma amplitude
% Voytek 
results.voytek = cfc_est_voytek(signal.theta_phase, signal.gamma_amp_theta);
