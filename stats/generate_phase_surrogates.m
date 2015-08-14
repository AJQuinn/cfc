function surrogates = generate_phase_surrogates( data, n, mean_term, std_term )
%% Generate surrogate data though method 3 in Hurtado et al 2004. Statistical
% method for detection of phase locking episodes in neural oscillations. J
% Neurophysiol. 10.1152/jn.00853.2003.
%
% This scrambles the phase spectrum whilst preserving the amplitude spectrum
%
% data is [1 x nsamples]
% n is an integer denoting the number of surrogates

if nargin < 4 || isempty(std_term)
    % If the output std isn't defined, we'll use the std of the data
    std_term = mean(std(data,[],2));
end

if nargin < 3 || isempty(mean_term)
    % Set surrogate mean to zero
    mean_term = 0;
end

% Get amplitude spectrum
amp_spec = abs( fft( data, [], 2 ) );
nsamples = size(amp_spec,2);
if mod(nsamples,2)
    % We have an even number of components
    stop_sample = (nsamples+1)/2;
else
    % We have an odd number of components
    stop_sample = nsamples/2 + 1;
end
amp_spec = amp_spec(2:stop_sample);

% Generate phase noise
noise = rand(n,size(amp_spec,2)) .* (2*pi);

% Make new phase scrabled spectrum
rand_spec = bsxfun(@times, exp(1i*noise), amp_spec);

% flippit
rand_spec = [rand_spec,fliplr(conj(rand_spec))];

% Add mean term
rand_spec = [ones(n,1) * mean_term * nsamples, rand_spec];

% Create new time_course
surrogates = ifft(rand_spec,[],2);

% Sanity check
if ~isreal(surrogates)
    warning('Surrogate time series are not real! something has gone wrong')
end

% Normalise time_series
 surrogates = bsxfun(@rdivide,surrogates, std(surrogates,[],2).*std_term );

