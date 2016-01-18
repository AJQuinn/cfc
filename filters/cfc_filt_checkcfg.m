function filt_cfg = cfc_filt_checkcfg( filt_cfg )
%function filt_cfg = cfc_filt_checkcfg( filt_cfg )
%
% Sanity check filter parameters, raise warnings if something silly might
% happen. Add default values if they aren't provided

%% Set transband if not defined
if ~isfield(filt_cfg,'trans_width')
    warning('Using default transition bandwidth')
    filt_cfg.trans_width = filt_cfg.pass_width + 10;
end

%% Check we aren't too close to nyquist
nyquist = filt_cfg.sample_rate/2;
if isfield(filt_cfg,'centre_freq')

    if filt_cfg.centre_freq >= nyquist
        warning( 'Centre frequency is above Nyquist!')
    end
    if filt_cfg.centre_freq + (filt_cfg.pass_width/2) >= nyquist
       warning( 'Upper pass band is above Nyquist!')
    end
    if filt_cfg.centre_freq + (filt_cfg.trans_width/2) >= nyquist
       warning( 'Upper transition band is above Nyquist!')
    end
end

if isfield(filt_cfg,'explore_freq')

    if length(filt_cfg.explore_freq) == 2
        freq = filt_cfg.explore_freq(2);
    else
        freq = filt_cfg.centre_freq;
    end
    
    if freq >= nyquist
        warning( 'Upper explore frequency is above Nyquist!')
    end
    if freq + (filt_cfg.pass_width/2) >= nyquist
       warning( 'Upper pass band is above Nyquist!')
    end
    if freq + (filt_cfg.trans_width/2) >= nyquist
       warning( 'Upper transition band is above Nyquist!')
    end
end

%% Check we aren't too close to zero
if isfield(filt_cfg,'centre_freq')

    if filt_cfg.centre_freq <= 0
        warning( 'Centre frequency is below zero!')
    end
    if filt_cfg.centre_freq - (filt_cfg.pass_width/2) <= 0
       warning( 'Lower pass band is below zero!')
    end
    if filt_cfg.centre_freq - (filt_cfg.trans_width/2) <= 0
       warning( 'Lower transition band is below zero!')
    end
end

if isfield(filt_cfg,'explore_freq')

    if length(filt_cfg.explore_freq) == 2
        freq = filt_cfg.explore_freq(1);
    else
        freq = filt_cfg.centre_freq;
    end
    
    if freq <= 0
        warning( 'Lower explore frequency is below zero')
    end
    if freq - (filt_cfg.pass_width/2) <= 0
       warning( 'Lower pass band is below zero!')
    end
    if freq - (filt_cfg.trans_width/2) <= 0
       warning( 'Lower transition band is below zero!')
    end
end

%%
