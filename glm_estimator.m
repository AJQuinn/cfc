function [glm] = glm_estimator(modulating_signal,modulated_signal, distr)

% Function for estimating the GLM phase amplitude coupling as defined in Penny et al 2008
% Modulating signal is the theta phase and the modulated signal is the gamma amplitude

if nargin < 3
    distr = 'normal';
end

% Create design matrix
X = [ cos(modulating_signal); ...
      sin(modulating_signal); 
      ones(size(modulating_signal)) ]';
  
[beta,dev,stats] = glmfit(X,modulated_signal',distr,'constant','off');

glm = [];
glm.beta = beta;
glm.dev = dev;
glm.stats = stats;

ss_data = sum( modulated_signal .* modulated_signal);
ss_resid = sum ( stats.resid .* stats.resid );

glm.r2 = (ss_data - ss_resid) / ss_data;

