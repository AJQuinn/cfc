function [out] = cfc_peak_query(indata, cfg)
%% cfc_peak_detect
% Function using Savitsky-Golay smoothing filter to detect peaks in the
% spectrum of a given time series
%
% Inputs
% ------
% indata: cell array of 2d arrays
%       cell array containing data arrays to find peaks in [channels, samples],
%       will average over channels if more than one
% sample_rate: double
%       sampling frequency of the data
% freq_of_interest: double
%       frequency of interest to test across participants
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
        error('Please specify a frequency of interest in Hz');
    end

    if ~isfield(cfg, 'sample_rate')
        cfg.cfg.sample_rate = 1;
    end

    if ~isfield(cfg,'input_domain')
        error('Please specify the domain of the input data. Either time or frequency');
    end

    nppts = length(indata);
    out.ppt = cell(nppts,1);

    for ippt = 1:nppts

        obj = [];
        data = indata{ippt};

        %% Smoothing window size for the SG filter
        sg_win = round(length(data)/100);
        if rem(sg_win,2) == 0
            sg_win = sg_win+1;
        end

        if strcmp(cfg.input_domain, 'time')
            % We have time domain data, fft and smooth it

            % Smoothing window size for the SG filter
            sg_win = round(length(data)/100);
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

            obj = get_freq_of_interest(obj, cfg.freq_of_interest);

        else
            error('Please define input domain as EITHER time or frequency')
        end

        if ippt == 1;
            tmp_idx = obj.freq_vect > 0 & obj.freq_vect < freq_of_interest+10;
            smo_psd = obj.smo_psd(tmp_idx);
            freq_vect = obj.freq_vect(tmp_idx);
        else
            smo_psd = cat(1,smo_psd,obj.smo_psd(tmp_idx));
        end

        %% Extract frequencies of interest
        halfwidth = .5;
        freq_idx_of_interest = obj.freq_vect > freq_of_interest-halfwidth & obj.freq_vect < freq_of_interest+halfwidth;
        obj.fft = obj.fft(freq_idx_of_interest);
        obj.freq_vect = obj.freq_vect(freq_idx_of_interest);
        obj.smo_psd = obj.smo_psd(freq_idx_of_interest);
        if detrend == true
            obj.sub_vect = obj.sub_vect(freq_idx_of_interest);
        end
        [~,freq_to_test_idx] = min(abs(obj.freq_vect == freq_of_interest));

        %% Subtract 1/f if requested
        if strcmp(detrend,'1/f')
            obj.sub_vect = (1./obj.freq_vect) * obj.fft(1);
            obj.fft = obj.fft - obj.sub_vect;
            obj.smo_psd = obj.smo_psd - obj.sub_vect;
        elseif strcmp(detrend,'linear')
            X = obj.freq_vect;
            X = cat(1,X,ones(size(X)));
            Y_hat = X'* ( X'\obj.fft');
            obj.fft = obj.fft - Y_hat';
            Y_hat = X'* ( X'\obj.smo_psd');
            obj.smo_psd = obj.smo_psd - Y_hat';
        elseif strcmp(detrend,'polyfit')
            S.detrend_pred = zeros(size(obj.smo_psd,1),size(obj.smo_psd,2));
            modelfun = @(p,x)(p(1)*(x-p(4)) - p(2)*(x-p(4)) + p(3)*(x-p(4)).^2);
            beta0 = [1, 1, 1, 500];
            [beta,R,J,covb,mse] = nlinfit(obj.freq_vect,obj.smo_psd,modelfun,beta0);
            % Get confidence intervals for the parameters in the fit
            ci = nlparci(beta,R,'covar',covb);
            pred = modelfun(beta,obj.smo_psd);
            obj.smo_psd = obj.smo_psd - pred;
            %obj.smo_psd = pred;
            %S.data(:,idx) = S.data(:,idx) + abs(min(S.data(:,idx))); % no values less than 0
            figure;plot(pred)
            figure;plot(obj.smo_psd)
        end


        %% Perform a linear regression to the differential of the smoothed spectrum
        pk_cnt = 1;

        % Assign peak location and height
        tmp.loc = freq_of_interest;
        tmp.peak_freq = freq_of_interest;
        tmp.peak_amp = obj.smo_psd(freq_to_test_idx);

        % Create design matrix for regression
        %obj.freq_vect(obj.freq_vect > obj.freq_vect(lo)-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+fre
        tmp.X = obj.freq_vect;
        tmp.X = cat(1,tmp.X,ones(size(tmp.X)));

        % Take differential of the smoothed fft
        %smo_tmp = sgolayfilt(diff(obj.smo_psd),order,sg_win);
        smo_tmp = [0 diff(obj.smo_psd)];

        % Assign differential as predicted variable for regression
        %tmp.Y = smo_tmp(obj.freq_vect > obj.freq_vect(locs(ipk))-freq_width & obj.freq_vect < obj.freq_vect(locs(ipk))+freq_width);
        tmp.Y = smo_tmp;

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

        % Store results
        obj.peak{1} = tmp;
        out.ppt{ippt} = obj;

    end

    % Store some key information at the top level
    out.peak_amplitudes = arrayfun(@(i) out.ppt{i}.peak{1}.peak_amp, 1:nppts);
    out.peak_frequencies = arrayfun(@(i) out.ppt{i}.peak{1}.peak_freq, 1:nppts);
    [~,out.peaks_by_amplitude] = sort(out.peak_amplitudes,2,'descend');

    %% Group level test
    betas = arrayfun(@(i) out.ppt{i}.peak{1}.coeff(1), 1:nppts);
    R2 = arrayfun(@(i) out.ppt{i}.peak{1}.r2, 1:nppts);
    SE = arrayfun(@(i) out.ppt{i}.peak{1}.SE, 1:nppts);
    tvals = arrayfun(@(i) out.ppt{i}.peak{1}.t, 1:nppts);

    % currently fixed effects, would be awfully nice to use SE here...
    [~,P,ci,stats] = ttest(tvals);
    out.P = P;
    out.ci = ci;
    out.stats = stats;

    %% Make a fancy plot


    figure
    subplot(211)
    plot(freq_vect,smo_psd);hold on
    plot(freq_vect,mean(smo_psd,1),'linewidth',3)
    subplot(212)
    plot(freq_vect(2:end),diff(smo_psd'));hold on
    for idx = 1:nppts
        plot(out.ppt{idx}.peak{1}.X(1,:),out.ppt{idx}.peak{1}.Y_hat,'k');
    end

    %% Print a fancy table
    header = sprintf('%-8s%-10s%-12s%-12s%-8s%-8s%-8s%-8s','Peak','Freq(Hz)','Amp','Beta','R2','SE','T','p');
    table_lines = repmat('-',1,length(header));

    disp(table_lines);disp(header);disp(table_lines);
    for idx = 1:nppts
        fprintf('%-8d%-10.2f%-12.2f%-12.2f%-8.2f%-8.2f%-8.2f%-8.5f\n',...
            idx,...
            out.ppt{idx}.peak{1}.peak_freq,...
            out.ppt{idx}.peak{1}.peak_amp,...
            out.ppt{idx}.peak{1}.coeff(1),...
            out.ppt{idx}.peak{1}.r2,...
            out.ppt{idx}.peak{1}.SE,...
            out.ppt{idx}.peak{1}.t,...
            out.ppt{idx}.peak{1}.p);
    end
    disp(table_lines);
    fprintf('%-8s%-10.2f%-12.2f%-12.2f%-8.2f%-8.2f%-8.2f%-8.5f\n',...
        'total',...
        mean(out.peak_frequencies),...
        mean(out.peak_amplitudes),...
        mean(betas),...
        nan,nan,...
        stats.tstat,...
        P);
    disp(table_lines);


end
