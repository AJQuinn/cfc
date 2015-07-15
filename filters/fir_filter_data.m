function [filt_data,D] = fir_filter_data(data,filter_cfg,padding,method)
%
% filter_cfg: struct containing
%   order:
%   sample_rate:
%   centre_freq:
%   pass_width:
%   trans_width:

if nargin < 4
    % Onepass is the more computationally efficient options as, for a fir filter, the
    % phase delay can be corrected with a simple linear time-shift.
    % Two pass uses filtfilt which runs the filter twice to remove any
    % phase lag
    method = 'onepass';
end

if nargin < 3
    padding = round(filter_cfg.sample_rate);
end

%% Make filter
[D,~,~,~] = make_filter(filter_cfg);

%% Pad the data
data = [zeros(padding,1); data; zeros(padding,1)];

%% Filter the data

if isa(D,'digitalFilter')
    % Filter using the object
    
    if strcmp(method,'twopass')
        [filt_data] = filtfilt(D,data);
        filt_data = filt_data(padding:end-padding);
    elseif strcmp(method,'onepass')
        [filt_data] = filter(D,data);
        pd = round(phasedelay(D));
        filt_data = filt_data(padding+pd(2):(end-padding-1)+pd(2));
    end
    
else
    % Filter using the vector
    
    if strcmp(method,'twopass')
        [filt_data] = filtfilt(D,1,data);
        filt_data = filt_data(padding:end-padding);
    elseif strcmp(method,'onepass')
        [filt_data] = filter(D,1,data);
        pd = round(phasedelay(D,1));
        filt_data = filt_data(padding+pd(2):(end-padding-1)+pd(2));
    end
    
end

