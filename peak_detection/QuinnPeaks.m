function [obj] = QuinnPeaks(data, sample_rate, freq_of_interest, order,detrend,fft_len)
%% QuinnPeaks
% Function using Savitsky-Golay smoothing filter to detect peaks in the
% spectrum of a given time series
%
% Inputs
% ------
% data: 2d array
%       series to find peaks in [channels, samples], will average over channels
%       if more than one
% sample_rate: double
%       sampling frequency of the data
% freq_of_interest: vector
%       low and high frequencies of interest eg [ .01 100 ]
% order: int [optional]
%       order for sgolayfilter
% detrend: str [optional]
%       Optional detrending of spectrum. Choose from '1/f','linear' or 'polyfit'

    if nargin < 6 || isempty(fft_len)
        fft_len = [];
    end

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
    obj.fft = fftshift(abs(fft(data,fft_len,2)));
    obj.freq_vect = linspace(-sample_rate/2,sample_rate/2,length(obj.fft));

    if size(data,1) > 1
        obj.fft = mean(obj.fft,1);
    end

    %% Smooth spectrum
    disp('Smoothing')
    obj.smo_fft = sgolayfilt(obj.fft, order, sg_win);

    %% Extract frequencies of interest
    freq_idx_of_interest = obj.freq_vect > freq_of_interest(1) & obj.freq_vect < freq_of_interest(2);
    obj.fft = obj.fft(freq_idx_of_interest);
    obj.freq_vect = obj.freq_vect(freq_idx_of_interest);
    obj.smo_fft = obj.smo_fft(freq_idx_of_interest);
    if detrend == true
        obj.sub_vect = obj.sub_vect(freq_idx_of_interest);
    end

    %% Subtract 1/f if requested
    if strcmp(detrend,'1/f')
        obj.sub_vect = (1./obj.freq_vect) * obj.fft(1);
        obj.fft = obj.fft - obj.sub_vect;
        obj.smo_fft = obj.smo_fft - obj.sub_vect;
    elseif strcmp(detrend,'linear')
        X = obj.freq_vect;
        X = cat(1,X,ones(size(X)));
        Y_hat = X'* ( X'\obj.fft');
        obj.fft = obj.fft - Y_hat';
        Y_hat = X'* ( X'\obj.smo_fft');
        obj.smo_fft = obj.smo_fft - Y_hat';
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



    %% Fit diff of peaks
    disp('Peaking')
    freq_width = 1;
    % Find peaks using file-exchange script as matlab's one sucks
    % algorithm from: http://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder
    % Peaks must be larger than the median differetial step in the spectrum
    [locs,pks] = peakfinder(obj.smo_fft,median(diff(obj.smo_fft)));
    % Remove peaks within 5 samples of the start or end
    locs(locs < 5 | locs > size(obj.freq_vect,2)-5) = [];
    pks(locs < 5 | locs > size(obj.freq_vect,2)-5) = []

    %% Perform a linear regression to the differential of the smoothed spectrum
    % +/- 1Hz around each peak in the spectrum, this might be a bit big, .5Hz is probably fine...
    disp('Peak ID')
    pk_cnt = 1;
    for ipk = 1:length(pks)

        % Assign peak location and height
        tmp.loc = locs(ipk);
        tmp.peak_amp = pks(ipk);

        % Create design matrix for regression
        tmp.X = obj.freq_vect(obj.freq_vect > obj.freq_vect(locs(ipk))-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+freq_width);
        tmp.X = cat(1,tmp.X,ones(size(tmp.X)));

        % Take differential of the smoothed fft
        %smo_tmp = sgolayfilt(diff(obj.smo_fft),order,sg_win);
        smo_tmp = [0 diff(obj.smo_fft)];

        % Assign differential as predicted variable for regression
        tmp.Y = smo_tmp(obj.freq_vect > obj.freq_vect(locs(ipk))-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+freq_width);

        % Regress!
        tmp.coeff = tmp.X'\tmp.Y';

        % Get predicted spectrum
        tmp.Y_hat = tmp.X'*tmp.coeff;

        % Get standard error sqrt [ ?(yi - ?i)2 / (n - 2) ] / sqrt [ ?(xi - x)2 ]
        denom = sqrt ( sum ( (tmp.X(1,:) - mean(tmp.X(1,:))) .^2 ) );
        num   = sqrt ( sum ( (tmp.Y - tmp.Y_hat').^2 ) / ( size(tmp.Y,2) - 2));
        tmp.SE = num / denom;

        % Get t-value for peak (testing against a beta of 0)
        tmp.t = tmp.coeff(1) / tmp.SE;

        % Get p-value from t-distribution
        tmp.p = 2*(1-tcdf(abs(tmp.t),size(tmp.Y,2)-2));

        % Get R2
        tmp.r2 = 1 - ( sum((tmp.Y - tmp.Y_hat').^2) / sum(tmp.Y.^2) );

        % Create an interpolated oj for sub-sample resolution
        interp_X = linspace(tmp.X(1,1),tmp.X(1,end),10000);
        tmp.interp_X = cat(1,interp_X,ones(size(interp_X)));
        tmp.interp_Y_hat = tmp.interp_X'*tmp.coeff;

        % peak frequency is then the zero crossing
        zero_crossing = tmp.interp_X(1,diff(sign(tmp.interp_Y_hat)) ~= 0)
        if isempty(zero_crossing)
            tmp.peak_freq = nan;
        else
            tmp.peak_freq = tmp.interp_X(1,diff(sign(tmp.interp_Y_hat)) ~= 0);
            % Save in main structure
            obj.peak{pk_cnt} = tmp;
            pk_cnt = pk_cnt+1;
        end

    end

    pk_cnt = pk_cnt - 1
    
    % Store some key information at the top level
    obj.peak_amplitudes = arrayfun(@(i) obj.peak{i}.peak_amp, 1:pk_cnt);
    obj.peak_frequencies = arrayfun(@(i) obj.peak{i}.peak_freq, 1:pk_cnt);
    [~,obj.peaks_by_amplitude] = sort(obj.peak_amplitudes,2,'descend');

    %% Print a fancy table
    header = sprintf('%-8s%-10s%-12s%-12s%-8s%-8s%-8s%-8s','Peak','Freq(Hz)','Amp','Beta','R2','SE','T','p');
    table_lines = repmat('-',1,length(header));

    disp(table_lines);disp(header);disp(table_lines);
    for idx = 1:pk_cnt
        fprintf('%-8d%-10.2f%-12.2f%-12.2f%-8.2f%-8.2f%-8.2f%-8.2f\n',...
            idx,...
            obj.peak{idx}.peak_freq,...
            obj.peak{idx}.peak_amp,...
            obj.peak{idx}.coeff(1),...
            obj.peak{idx}.r2,...
            obj.peak{idx}.SE,...
            obj.peak{idx}.t,...
            obj.peak{idx}.p);
    end
    disp(table_lines);

end
