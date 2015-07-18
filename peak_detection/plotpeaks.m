function plot_sgolay_peaks(obj, outpath, max_freq)
%%
% Create a plot from the output struct from sgolay_peaks
% 

if nargin<3
    max_freq = obj.freq_vect(end-1);
end

if nargin<2
    outpath =false
end

max_freq_idx = length(obj.freq_vect(obj.freq_vect < max_freq));

figure;
subplot(211)
plot(obj.freq_vect(1:max_freq_idx),obj.fft(1:max_freq_idx))
hold on
plot(obj.freq_vect(1:max_freq_idx),obj.smo_fft(1:max_freq_idx),'r')

for ipk = 1:length(obj.peak)
    plot(obj.freq_vect(obj.peak{ipk}.loc),obj.peak{ipk}.peak,'or');
end

grid on;


%% Smoothing window size for the SG filter
sg_win = round(length(obj.freq_vect)/100);
if rem(sg_win,2) == 0
    sg_win = sg_win+1;
end
smo_tmp = sgolayfilt(diff(obj.smo_fft),2,sg_win);


subplot(212)
plot(obj.freq_vect(2:max_freq_idx),diff(obj.smo_fft(1:max_freq_idx)))
hold on
plot(obj.freq_vect(2:max_freq_idx),smo_tmp(1:max_freq_idx-1),'r')
grid on;

for ipk = 1:length(obj.peak)
    plot(obj.peak{ipk}.X(1,:),obj.peak{ipk}.diff_obj,'g','linewidth',2)
end

if outpath ~= false
    savefig(outpath);
    print(gcf, outpath, '-deps');
end

