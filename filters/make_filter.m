function [D,h,phi,w] = make_filter(varargin)
%
% USAGE 1:
% Call either with a single argument with a filter_cfg object
%
% [D,h,phi,w] = make_filter(filter_cfg)
%
% in which filter_cfg is a struct containing
%   order:
%   sample_rate:
%   centre_freq:
%   pass_width:
%   trans_width:
%
% USAGE 2:
% [D,h,phi,w] = make_filter(order, sample_rate, ...
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
    warning('Wrong number of arguments passed, check make_filter usage');
end
        
if centre_freq-(trans_width/2) <= 0
    warning('Lower filter transition band is too close to zero!');
end

if centre_freq+(trans_width/2) >= sample_rate/2
    warning('Upper filter transition band is too close to Nyquist!')
end

%% Design a filter
D = designfilt('bandpassfir', 'FilterOrder', order, ...
      'StopbandFrequency1', centre_freq-(trans_width/2),...
      'PassbandFrequency1', centre_freq-(pass_width/2) , ...
      'PassbandFrequency2', centre_freq+(pass_width/2),...
      'StopbandFrequency2', centre_freq+(trans_width/2),...
      'SampleRate', sample_rate);
  
  % Get frequency response
  [h,~] = freqz(D);
  
  % Get phase response
  [phi,w] = phasez(D);