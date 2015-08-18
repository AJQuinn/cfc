function [mean_vector] = mi_norm_estimator(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006. The MI
% is returned as a z-value on a null distribution created by phase shifted
% surrogates. See supplimentary materials of Canolty et al 2006.

[~,nsamples] = size(modulating_signal);

% Combine into complex form
Z = modulated_signal .* exp(1i * modulating_signal);

% Compute mean vector
mean_vector = mean(squeeze(Z),2);

% Compute surrogates
nsurrogates = 200;
surrogate_Z = nan(nsurrogates,1);
skips = randi([fix(nsamples/10),9*fix(nsamples/10)],nsurrogates,1);
for idx = 1:nsurrogates
    surrogate_Z(idx) = abs(mean( circshift(modulated_signal,skips(idx),2).*exp(1i * modulating_signal) ));
end

% Normalise
Z_length = (abs(mean_vector)-mean(surrogate_Z))./std(surrogate_Z);
mean_vector = abs(Z_length.*exp(1i*angle(mean_vector)));