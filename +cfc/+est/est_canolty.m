function [m_norm] = est_canolty(modulating_signal,modulated_signal,srate)

penny = true; % Whether we want M_norm or M_raw, penny suggests M_raw

numpoints=length(modulating_signal); %% number of sample points in raw signal

if penny == false
    numsurrogate=200; %% number of surrogate values to compare to actual value
    minskip=srate; %% time lag must be at least this big
    maxskip=numpoints-srate; %% time lag must be smaller than this
    skip=ceil(numpoints.*rand(numsurrogate*2,1));
    skip(find(skip>maxskip))=[];
    skip(find(skip<minskip))=[];
    skip=skip(1:numsurrogate,1);
    surrogate_m=zeros(numsurrogate,1);
end

amplitude = modulated_signal;
phase = modulating_signal;


%% complex-valued composite signal
z=amplitude.*exp(1i*phase);
%% mean of z over time, prenormalized value
m_raw=mean(z);

if penny == false
    %% compute surrogate values
    for s=1:numsurrogate
        surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)];
        surrogate_m(s)=abs(mean(surrogate_amplitude.*exp(1i*phase)));
        %disp(numsurrogate-s)
    end
    %% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
    [surrogate_mean,surrogate_std]=normfit(surrogate_m);
    %% normalize length using surrogate data (z-score)
    m_norm_length=(abs(m_raw)-surrogate_mean)/surrogate_std;
    
    m_norm_phase=angle(m_raw);
    m_norm=m_norm_length*exp(1i*m_norm_phase);
   
else
    m_norm = m_raw;
end