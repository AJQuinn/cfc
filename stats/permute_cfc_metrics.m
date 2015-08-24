function [results] = permute_cfc_metrics(signal,alpha,nperms)
% Window permutations to establsh a null distribution

%% Housekeeping
results = struct('plv_null',    zeros(nperms,1),...
    'aec_null',    zeros(nperms,1),...
    'esc_null',    zeros(nperms,1),...
    'nesc_null',   zeros(nperms,1),...
    'glm_null',    zeros(nperms,1),...
    'mi_null',     zeros(nperms,1));

msg = '';
%% Estimate null distributions
for idx = 1:nperms
    
    fprintf(repmat(char(8),1,length(msg)));
    msg = sprintf('Calculating permutation %d/%d', ...
        idx,nperms);
    fprintf(msg);
    
    signal = permute_signal(signal,'window');
    
    % Estimate the PLV - using theta phase and gamma amplitude
    % Lachaux 1999
    results.plv_null(idx) = max(plv_estimator(signal.theta_phase, signal.gamma_amp_phase));
    
    % Estimate AEC - using theta amplitude and gamma amplitude
    % Bruns & Eckhorn 2004
    results.aec_null(idx) = max(cfc_est_aec(signal.theta_amp, signal.gamma_amp));
    
    % Estimate MI - using theta phase and gamma amplitude
    % Canolty 2006
    results.mi_null(idx) = max(mi_estimator(signal.theta_phase,signal.gamma_amp));
    
    % Estimate ESC - using theta signal and gamma amplitude
    % Bruns & Eckhorn 2004
    results.esc_null(idx) = max(esc_estimator(signal.theta,signal.gamma_amp));
    
    % Estimate NESC - using theta phase and gamma amplitude
    % Penny 2007
    results.nesc_null(idx) = max(nesc_estimator(signal.theta_phase, signal.gamma_amp));
    
    % Estimate GLM -  using theta phase and gamma amplitude
    % Penny 2007
    results.glm_null(idx) = max(glm_estimator(signal.theta_phase, signal.gamma_amp));
    
end

%% Generate critcal values
results.alpha = alpha;

results.plv_cv = prctile(results.plv_null,1-alpha);
results.aec_cv = prctile(results.aec_null,1-alpha);
results.mi_cv = prctile(results.mi_null,1-alpha);
results.esc_cv = prctile(results.esc_null,1-alpha);
results.nesc_cv = prctile(results.nesc_null,1-alpha);
results.glm_cv = prctile(results.glm_null,1-alpha);

