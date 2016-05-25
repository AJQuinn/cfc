function [D,h,phi,w] = cfc_filt_generate(varargin)
%
% USAGE 1:
% Call either with a single argument with a filter_cfg object
%
% [D,h,phi,w] = cfc_filt_generate(filter_cfg)
%
% in which filter_cfg is a struct containing
%   order:
%   sample_rate:
%   centre_freq:
%   pass_width:
%   trans_width:
%
% USAGE 2:
% [D,h,phi,w] = cfc_filt_generate(order, sample_rate, ...
%                           centre_freq, pass_width, trans_width)

%% Housekeeping
if length(varargin) == 1
    % assume first argument is a filter_cfg
    obj = varargin{1};
    order = obj.order;
    sample_rate = obj.sample_rate;
    centre_freq = obj.centre_freq;
    pass_width = obj.pass_width;
    trans_width = obj.trans_width;
elseif length(varargin) == 5
    order = varargin{1};
    sample_rate = varargin{2};
    centre_freq = varargin{3};
    pass_width = varargin{4};
    trans_width = varargin{5};
else
    warning('Wrong number of arguments passed, check cfc_filt_generate usage');
end
        

pass = [centre_freq-(pass_width/2) centre_freq+(pass_width/2)];
trans = [centre_freq-(trans_width/2) centre_freq+(trans_width/2)];



if trans(1) <= 0
    %warning('Lower filter transition band is too close to zero!');
    trans(1) = 0;
end

if trans(2) >= sample_rate/2
    warning('Upper filter transition band is too close to Nyquist!')
    trans(2) = sample_rate/2;
end

% TODO: expose this option?
method = 'parks';

if strcmp(method, 'window')
    %% Design a filter with windowing method
    D = designfilt('bandpassfir', 'FilterOrder', order, ...
        'StopbandFrequency1', trans(1),...
        'PassbandFrequency1', pass(1) , ...
        'PassbandFrequency2', pass(2),...
        'StopbandFrequency2', trans(2),...
        'SampleRate', sample_rate);
    
    % Get frequency response
    [h,~] = freqz(D);
    
    % Get phase response
    [phi,w] = phasez(D);
    
else
    %% Design a filter with Parks & McClelland
    
    f = [0, trans(1)....
        pass(1),...
        pass(2),...
        trans(2),...
        sample_rate / 2];
    f = f / ( sample_rate / 2); % normalised frequency
    a = [0 0.0 1.0 1.0 0.0 0];
    D = firpm(order,f,a);
    
    % get frequency response
    [h,~] = freqz(D,1,512);
    
    % get phase response
    [phi,w] = phasez(D,1,512);
end
