function [C] = corr(modulating_signal,modulated_signal)
%function [C] = cfc.est.corr(modulating_signal,modulated_signal)
%
% Function for estimating the column-by-column correlation between two signals.
% This can be used to estimate several metrics when used with different inputs.
%
%   aec  = cfc_est_corr( theta_amp, gamma_amp )
%   esc  = cfc_est_corr( theta, gamma_amp )
%   nesc = cfc_est_corr( cos(theta), gamma_amp )
%
% Modulating signal is [1 x samples x realisations]
% Modulated signal is [1 x samples x realisations]
%
% A realisation can be either an epoch or a sliding window

if ndims(modulating_signal) == 2
    modulating_signal = modulating_signal';
    modulated_signal = modulated_signal';
else
    modulating_signal = squeeze(modulating_signal);
    modulated_signal = squeeze(modulated_signal);
end

%% Compute correlation
% http://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab

% Subtract mean column-wise
modulating_signal = bsxfun(@minus,modulating_signal,mean(modulating_signal,1));
modulated_signal  = bsxfun(@minus,modulated_signal,mean(modulated_signal,1));

% Normalise column-wise
modulating_signal = bsxfun(@times,modulating_signal,1./sqrt(sum(modulating_signal.^2,1)));
modulated_signal  = bsxfun(@times,modulated_signal,1./sqrt(sum(modulated_signal.^2,1)));

% Compute correlation
C=sum(modulating_signal.*modulated_signal,1);
