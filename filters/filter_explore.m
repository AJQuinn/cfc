function [] = filter_explore(obj)
%
%
% Input object
%
% pass_width
% trans_width
% lo_bound
% hi_bound
% freq_step
% sample_rate
% filter_order
% method: 'eeglab','designfilt'

S.obj = obj;

% Create figure
S.fh = figure('units','pixels',...        
              'position',[300 300 800 600],...
              'menubar','none',...
              'name','slider_plot',...
              'numbertitle','off',...
              'resize','off');  

% Design intial filter
if strcmp(obj.method,'designfilt')
  [~,h,phi,w] = make_filter(248,obj.sample_rate,obj.lo_bound,...
                            obj.pass_width,obj.trans_width);
elseif strcmp(obj.method,'eeglab')
  [~,filtwts] = eegfilt_silent(randn(1,1024),obj.sample_rate,...
                               obj.lo_bound-(obj.pass_width/2),...
                               obj.lo_bound+(obj.pass_width/2),...
                               0,[],0,'fir1');
  [h,w] = freqz(filtwts);
  [phi,w] = phasez(filtwts);
elseif strcmp(obj.method,'ft')
    N = 3*fix(obj.sample_rate / obj.lo_bound);
   [B, A] = fir1(N, [obj.lo_bound/(obj.sample_rate/2) obj.hi_bound/(obj.sample_rate/2)]);
   [h,w] = freqz(B,1,512);
end
   
% Get frequency vector w = (2*pi*f) / fw
freq_vect = (w*obj.sample_rate) / (2*pi);


% Create main axis
S.ax1 = axes('unit','pix',...
            'position',[40 375 720 200]);hold on;
S.ax2 = axes('unit','pix',...
            'position',[40 150 720 200]);hold on;       
% Plot magnitude response
S.mag_resp = plot(S.ax1,freq_vect,20*log(abs(h))); 
S.phase_resp = plot(S.ax2,freq_vect,phi);

% Plot pass band
ylim = get(S.ax1,'ylim');
pass_resp = zeros(size(freq_vect,1),size(freq_vect,2)) + ylim(1);
pass_resp(freq_vect > obj.lo_bound-(obj.pass_width/2) & ...
          freq_vect < obj.lo_bound+(obj.pass_width/2)) = 0;
S.pass_resp = plot(S.ax1,freq_vect,pass_resp,'r--');

% Plot transition band
trans_resp = zeros(size(freq_vect,1),size(freq_vect,2)) + ylim(1);
trans_resp(freq_vect > obj.lo_bound-(obj.trans_width/2) & ...
           freq_vect < obj.lo_bound+(obj.trans_width/2)) = 0;
S.trans_resp = plot(S.ax1,freq_vect,trans_resp,'g--');

% Housekeeping
grid on;          
ylabel(S.ax1,'Magnitude (dB)');
title(S.ax1,'Magnitude Response');
title(S.ax2,'Phase Response');

xlabel(S.ax2,'Frequency (Hz)');

legend(S.ax1,{'Magnitude Response','Pass Band','Transition Band'});

S.text1=uicontrol('style','text',...
         'position',[40 100 360 20],...
         'String',...
         ['Passband:',num2str(obj.lo_bound-(obj.pass_width/2)),...
         ' ', num2str(obj.lo_bound+(obj.pass_width/2)) ]);
S.text2=uicontrol('style','text',...
         'position',[40 80 360 20],...
         'String',...
         ['Trans band:',num2str(obj.lo_bound-(obj.trans_width/2)),...
         ' ', num2str(obj.lo_bound+(obj.trans_width/2)) ]);
     
% Create slider
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[40 50 720 10],...
                 'min',obj.lo_bound,'max',obj.hi_bound,'val',obj.lo_bound,...
                 'sliderstep',[1/(obj.hi_bound-obj.lo_bound) 1/(obj.hi_bound-obj.lo_bound)],...
                 'callback',{@sl_call,S}); 
             
             
function [] = sl_call(varargin)
% Callback for the slider.
[H,S] = varargin{[1,3]};  % calling handle and data structure.

value=round(get(H,'Value'));
set(H,'Value',value);


% Design new filter
if strcmp(S.obj.method,'designfilt')
[~,h,phi,w] = make_filter(248,S.obj.sample_rate,get(H,'value'),...
                            S.obj.pass_width,S.obj.trans_width);
elseif strcmp(S.obj.method,'eeglab')
  [~,filtwts] = eegfilt_silent(randn(1,1024),S.obj.sample_rate, ...
                               get(H,'value')-(S.obj.pass_width/2),...
                               get(H,'value')+(S.obj.pass_width/2),...
                               0,[],0,'fir1');
  [h,w] = freqz(filtwts);
  [phi,w] = phasez(filtwts);
end
freq_vect = (w*S.obj.sample_rate) / (2*pi);

ylim = get(S.ax1,'ylim');
pass_resp = zeros(size(freq_vect,1),size(freq_vect,2)) + ylim(1);
pass_resp(freq_vect > get(H,'value')-(S.obj.pass_width/2) & ...
          freq_vect < get(H,'value')+(S.obj.pass_width/2)) = 0;

trans_resp = zeros(size(freq_vect,1),size(freq_vect,2)) + ylim(1);
trans_resp(freq_vect > get(H,'value')-(S.obj.trans_width/2) & ...
           freq_vect < get(H,'value')+(S.obj.trans_width/2)) = 0;

% Update the plots
set(S.mag_resp,'ydata',20*log(abs(h)));
set(S.phase_resp,'ydata',phi);

set(S.pass_resp,'ydata',pass_resp);
set(S.trans_resp,'ydata',trans_resp);

set(S.text1,'String',['Passband:',num2str(get(H,'value')-(S.obj.pass_width/2)),...
         ' ', num2str(get(H,'value')+(S.obj.pass_width/2)) ]);
set(S.text2,'String',['Transband:',num2str(get(H,'value')-(S.obj.trans_width/2)),...
         ' ', num2str(get(H,'value')+(S.obj.trans_width/2)) ]);