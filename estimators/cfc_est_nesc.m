function [nesc] = ncfc_est_esc(modulating_signal,modulated_signal)

% Function for estimating the normalised Envelope-to-signal correlation as defined in Penny et al 2008
% Modulating signal is the theta filtere time series and modulated_signal is the gamma amplitude

modulating_signal = cos(modulating_signal)';
modulated_signal = modulated_signal';

%% Compute correlation
% http://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab

% Subtract mean column-wise
modulating_signal = bsxfun(@minus,modulating_signal,mean(modulating_signal,1));
modulated_signal  = bsxfun(@minus,modulated_signal,mean(modulated_signal,1));

% Normalise column-wise
modulating_signal = bsxfun(@times,modulating_signal,1./sqrt(sum(modulating_signal.^2,1)));
modulated_signal  = bsxfun(@times,modulated_signal,1./sqrt(sum(modulated_signal.^2,1)));

% Compute correlation
nesc=sum(modulating_signal.*modulated_signal,1);
