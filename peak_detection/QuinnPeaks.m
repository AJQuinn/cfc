function [obj] = QuinnPeaks(data, sample_rate, freq_of_interest, order,detrend)
%% Function using Savitsky-Golay smoothing filter to detect peaks in a
% spectrum or time-series
%
% data: vector
%       series to find peaks in
% sample_rate: scalar
%       sampling frequency of the data
% freq_of_interest: vector
%       low and high frequencies of interest eg [ .01 100 ]
% plot_flag: bool
%       flag to indicate whether to produce a plot


    if nargin < 5 || isempty(detrend)
        detrend = 'nope';
    end

    if nargin < 4 || isempty(order)
        order = 3;
    end

    if nargin < 3 || isempty(freq_of_interest)
        freq_of_interest = [ 0.001 100 ];
    end

    if nargin < 2 || isempty(sample_rate)
        sample_rate = 1;
    end

    obj = [];

    %% Smoothing window size for the SG filter
    sg_win = round(length(data)/100);
    if rem(sg_win,2) == 0
        sg_win = sg_win+1;
    end

    %% Peform frequency transform
    disp('FFT')
    obj.fft = fftshift(abs(fft(data,[],2)));
    obj.freq_vect = linspace(-sample_rate/2,sample_rate/2,length(obj.fft));

    
    if size(data,1) > 1
        obj.fft = mean(obj.fft,1);
    end

    %% Smooth spectrum
    disp('Smoothing')
    obj.smo_fft = sgolayfilt(obj.fft, order, sg_win);

    %% Subtract 1/f if requested
    if strcmp(detrend,'1/f')
        obj.sub_vect = (1./obj.freq_vect) * obj.fft(1);
        obj.fft = obj.fft - obj.sub_vect;
        obj.smo_fft = obj.smo_fft - obj.sub_vect;
    elseif strcmp(detrend,'polyfit')
        
        S.detrend_pred = zeros(size(obj.smo_fft,1),size(obj.smo_fft,2));
        modelfun = @(p,x)(p(1)*(x-p(4)) - p(2)*(x-p(4)) + p(3)*(x-p(4)).^2);
        beta0 = [1, 1, 1, 500];
        [beta,R,J,covb,mse] = nlinfit(obj.freq_vect,obj.smo_fft,modelfun,beta0);
        % Get confidence intervals for the parameters in the fit
        CI = nlparci(beta,R,'covar',covb);
        pred = modelfun(beta,obj.smo_fft);
        obj.smo_fft = obj.smo_fft - pred;
        %obj.smo_fft = pred;
        %S.data(:,idx) = S.data(:,idx) + abs(min(S.data(:,idx))); % no values less than 0
        figure;plot(pred)
        figure;plot(obj.smo_fft)
    end
        
        
    %% Extract frequencies of interest
    freq_idx_of_interest = obj.freq_vect > freq_of_interest(1) & obj.freq_vect < freq_of_interest(2);
    obj.fft = obj.fft(freq_idx_of_interest);
    obj.freq_vect = obj.freq_vect(freq_idx_of_interest);
    obj.smo_fft = obj.smo_fft(freq_idx_of_interest);
    if detrend == true
        obj.sub_vect = obj.sub_vect(freq_idx_of_interest);
    end
        
    %% Fit diff of peaks
    disp('Peaking')
    freq_width = .3;
    % Find peaks using file-exchange script as matlab's one sucks
    % algorithm from: http://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder
    %[locs,pks] = peakfinder(obj.smo_fft,max(obj.smo_fft)/50,[],[],false);
    [locs,pks] = peakfinder(obj.smo_fft,max(obj.smo_fft)/50,[],[],false);
    
    %% Perform a linear regression to the differential of the smoothed spectrum 
    % +/- 1Hz around each peak in the spectrum, this might be a bit big, .5Hz is probably fine...
    disp('Peak ID')
    for ipk = 1:length(pks)
        disp(ipk);
   
        % Assign peak location and height
        obj.peak{ipk}.loc = locs(ipk)
        obj.peak{ipk}.peak = pks(ipk)
 
        % Create design matrix for regression
        obj.peak{ipk}.X = obj.freq_vect(obj.freq_vect > obj.freq_vect(locs(ipk))-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+freq_width);
        obj.peak{ipk}.X = cat(1,obj.peak{ipk}.X,ones(size(obj.peak{ipk}.X)));

        % Take differential of the smoothed fft
        %smo_tmp = sgolayfilt(diff(obj.smo_fft),order,sg_win);
        smo_tmp = diff(obj.smo_fft);

        % Assign differential as predicted variable for regression
        obj.peak{ipk}.Y = smo_tmp(obj.freq_vect > obj.freq_vect(locs(ipk))-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+freq_width);
    
        % Regress!
        obj.peak{ipk}.coeff = obj.peak{ipk}.X'\obj.peak{ipk}.Y';
    
        % Get predicted spectrum
        obj.peak{ipk}.diff_obj = obj.peak{ipk}.X'*obj.peak{ipk}.coeff;

        %% Create an interpolated obj for sub-sample resolution
        interp_X = linspace(obj.peak{ipk}.X(1,1),obj.peak{ipk}.X(end,1),10000);
        obj.peak{ipk}.interp_X = cat(1,interp_X,ones(size(interp_X)));
        obj.peak{ipk}.interp_obj = obj.peak{ipk}.interp_X'*obj.peak{ipk}.coeff;

        % The zero crossing of this obj is the exact peak on x
        obj.peak{ipk}.zero_crossing = obj.peak{ipk}.interp_X(1,diff(sign(obj.peak{ipk}.interp_obj)) ~= 0);

        obj.peak{ipk}.interp_peak_loc = obj.freq_vect(locs(ipk));
        
    end

end