function [obj] = cfc_peak_detect( data, cfg )
%% cfc_peak_detect
% Function using Savitsky-Golay smoothing filter to detect peaks in the
% spectrum of a given time series
%
% Inputs
% ------
% data: 2d array
%       series to find peaks in [channels, samples], will average over channels
%       if more than one
%
% cfg struct defining options
%
% cfg.sample_rate: double
%       sampling frequency of the data
% cfg.freq_of_interest: vector
%       low and high frequencies of interest eg [ .01 100 ]
% cfg.input_domain: string
%       The domain of the input data, either 'time' or 'frequency'
% order: int [optional]
%       order for sgolayfilter
% detrend: str [optional]
%       Optional detrending of spectrum. Choose from '1/f','linear' or 'polyfit'

    if ~isfield(cfg, 'fft_len')
        cfg.fft_len = [];
    end

    if ~isfield(cfg, 'detrend')
        cfg.detrend = 'nope';
    end

    if ~isfield(cfg, 'smoothing_order')
        cfg.smoothing_order = 3;
    end

    if ~isfield(cfg, 'freq_of_interest')
        cfg.cfg.freq_of_interest = [ 0.001 100 ];
    end

    if ~isfield(cfg, 'sample_rate')
        cfg.cfg.sample_rate = 1;
    end

    if ~isfield(cfg,'input_domain')
        error('Please specify the domain of the input data. Either time or frequency');
    end

    obj = [];

    if strcmp(cfg.input_domain, 'time')
        % We have time domain data, fft and smooth it

        % Smoothing window size for the SG filter
        sg_win = round(length(data)/250);
        if rem(sg_win,2) == 0
            sg_win = sg_win+1;
        end

        % Peform frequency transform
        obj.psd = fftshift(abs(fft(data,cfg.fft_len,2)).^2);
        obj.freq_vect = linspace(-cfg.sample_rate/2,cfg.sample_rate/2,size(obj.psd,2));

        if size(data,1) > 1
            obj.psd = mean(obj.psd,1);
        end

        obj = get_freq_of_interest(obj, cfg.freq_of_interest);

        % Smooth spectrum
        obj.smo_psd = sgolayfilt(obj.psd, cfg.smoothing_order, sg_win);

    elseif strcmp(cfg.input_domain,'frequency')
        % We have frequency domain data, extract the interesting part

        obj.psd = data;
        obj.freq_vect = cfg.freq_vect;
        obj = get_freq_of_interest(obj, cfg.freq_of_interest);

        obj.smo_psd = obj.psd;

    else
        error('Please define input domain as EITHER time or frequency')
    end


    %% Subtract 1/f if requested
    if strcmp(cfg.detrend,'1/f')
        obj.sub_vect = (1./obj.freq_vect) * obj.psd(1);
        obj.psd = obj.psd - obj.sub_vect;
        obj.smo_psd = obj.smo_psd - obj.sub_vect;
    elseif strcmp(cfg.detrend,'linear')
        X = obj.freq_vect;
        X = cat(1,X,ones(size(X)));
        Y_hat = X'* ( X'\obj.psd');
        obj.psd = obj.psd - Y_hat';
        Y_hat = X'* ( X'\obj.smo_psd');
        obj.smo_psd = obj.smo_psd - Y_hat';
    elseif strcmp(cfg.detrend,'polyfit')
        S.detrend_pred = zeros(size(obj.smo_psd,1),size(obj.smo_psd,2));
        modelfun = @(p,x)(p(1)*(x-p(4)) - p(2)*(x-p(4)) + p(3)*(x-p(4)).^2);
        beta0 = [1, 1, 1, 500];
        [beta,R,J,covb,mse] = nlinfit(obj.freq_vect,obj.smo_psd,modelfun,beta0);
        % Get confidence intervals for the parameters in the fit
        CI = nlparci(beta,R,'covar',covb);
        pred = modelfun(beta,obj.smo_psd);
        obj.smo_psd = obj.smo_psd - pred;
        %obj.smo_psd = pred;
        %S.data(:,idx) = S.data(:,idx) + abs(min(S.data(:,idx))); % no values less than 0
        figure;plot(pred)
        figure;plot(obj.smo_psd)
    end

    %% Fit diff of peaks
    disp('Peaking')
    freq_width = 1;
    % Find peaks using file-exchange script as matlab's one sucks
    % algorithm from: http://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder
    % Peaks must be larger than the median differetial step in the spectrum
    [locs,pks] = peakfinder(obj.smo_psd,median(diff(obj.smo_psd)));
    % Remove peaks within 5 samples of the start or end
    pks(locs < 5 | locs > size(obj.freq_vect,2)-5) = [];
    locs(locs < 5 | locs > size(obj.freq_vect,2)-5) = [];

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
        %smo_tmp = sgolayfilt(diff(obj.smo_psd),order,sg_win);
        smo_tmp = [0 diff(obj.smo_psd)];

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
        zero_crossing = tmp.interp_X(1,diff(sign(tmp.interp_Y_hat)) ~= 0);
        if isempty(zero_crossing)
            tmp.peak_freq = nan;
        else
            tmp.peak_freq = tmp.interp_X(1,diff(sign(tmp.interp_Y_hat)) ~= 0);
            % Save in main structure
            obj.peak{pk_cnt} = tmp;
            pk_cnt = pk_cnt+1;
        end

    end

    pk_cnt = pk_cnt - 1;

    % Store some key information at the top level
    obj.peak_amplitudes = arrayfun(@(i) obj.peak{i}.peak_amp, 1:pk_cnt);
    obj.peak_frequencies = arrayfun(@(i) obj.peak{i}.peak_freq, 1:pk_cnt);

    % Return at least three peaks, use nans if there aren't enough.
    % This is to help the group plotting later...
    if length(obj.peak_amplitudes) < 3
        obj.peak_amplitudes = cat(2,obj.peak_amplitudes,ones(1,3-length(obj.peak_amplitudes))*nan);
        obj.peak_frequencies = cat(2,obj.peak_frequencies,ones(1,3-length(obj.peak_frequencies))*nan);
    end
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

function this_obj = get_freq_of_interest(this_obj, freq_of_interest)

    %% Extract frequencies of interest
    freq_idx_of_interest = this_obj.freq_vect > freq_of_interest(1) & this_obj.freq_vect < freq_of_interest(2);
    this_obj.psd = this_obj.psd(freq_idx_of_interest);
    this_obj.freq_vect = this_obj.freq_vect(freq_idx_of_interest);

end

