function [results] = estimate_cfc_metrics(signal)


results = [];

%% Estimate the PLV - using theta phase and gamma amplitude
% Lachaux 1999
results.plv = plv_estimator(signal.theta_phase, signal.gamma_amp_phase);

%% Estimate AEC - using theta amplitude and gamma amplitude 
% Bruns & Eckhorn 2004
results.aec= aec_estimator(signal.theta_amp, signal.gamma_amp);

%% Estimate MI - using theta phase and gamma amplitude
% Canolty 2006
results.mi = mi_estimator(signal.theta_phase,signal.gamma_amp);

%% Estimate ESC - using theta signal and gamma amplitude
% Bruns & Eckhorn 2004
results.esc = esc_estimator(signal.theta,signal.gamma_amp);

%% Estimate NESC - using theta phase and gamma amplitude
% Penny 2007
results.nesc = nesc_estimator(signal.theta_phase, signal.gamma_amp);

%% Estimate GLM -  using theta phase and gamma amplitude
% Penny 2007
results.glm = glm_estimator(signal.theta_phase, signal.gamma_amp);

%% Estimate Voytek - using theta phase and theta band filtered gamma amplitude
% Voytek 
results.voytek = voytek_estimator(signal.theta_phase, signal.gamma_amp_theta);
