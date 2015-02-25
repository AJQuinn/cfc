function [r2] = glm_estimator(modulating_signal,modulated_signal, distr)

% Function for estimating the GLM phase amplitude coupling as defined in Penny et al 2008
% Modulating signal is the theta phase and the modulated signal is the gamma amplitude

if nargin < 3
    distr = 'normal';
end

r2 = zeros(1,size(modulating_signal,2));
for ep = 1:size(modulating_signal,2)

    X = cat(2,cos(modulating_signal(:,ep)), ...
              sin(modulating_signal(:,ep)), ...
              ones(size(modulating_signal,1),1));
     
   [beta,dev,stats] = glmfit(X,modulated_signal(:,ep)',distr,'constant','off');

   glm = [];
   glm.beta = beta;
   glm.dev = dev;
   glm.stats = stats;

   ss_data = sum( modulated_signal(:,ep) .* modulated_signal(:,ep));
   ss_resid = sum ( stats.resid .* stats.resid );

   r2(ep) = (ss_data - ss_resid) / ss_data;

end

