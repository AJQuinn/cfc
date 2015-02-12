function [ data ] = make_sw_data(X,window_size,step)
%
% Function for splitting a [ samples x epochs ] dataset in to [ samples x windows x epochs ] separate sliding windows based on a window length and step size (must be set in samples rather than ms)

[nsamples,nepochs] = size(X);
nwindows = floor ( (nsamples - window_size) / step);
data = zeros(window_size,nwindows,nepochs);

for idx = 1:nwindows
    start_idx = (idx-1)*step + 1;
    end_idx = (idx-1)*step + window_size; 

    data(:,idx,:) = X(start_idx:end_idx,:);
end
