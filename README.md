# cfc
Cross frequency coupling toolbox in MatLab


# Getting started

Add the top directly to your matlab path, note that you do not need to use genpath or add any subdirectories by hand
```matlab
addpath('/Users/andrew/Software/cfc/')
```

You can then access the functions within the toolbox within a cfc struct. For instance
```matlab
cfc.example         % A set of example analyses using simulated data
cfc.simulate        % A function to generate simulated data
cfc.estimate_cfc    % The main function for estimating cross frequency coupling
```

Many analyses will require you to define a cfg struct to specify the parameters of the analysis. For instance to estimate Canolty's Modulation Index and the Phase Locking Value between 10 and 70Hz in a time series and compute 250 null permutations, you would define:
```matlab
cfc_cfg = struct('sr',      sample_rate,...
    'lo_freq',              10,...
    'lo_bandwidth',         2,...
    'hi_freq',              70,...
    'hi_bandwidth',         20,...
    'metrics',               {{'MI','PLV'}},...
    'nperms',               250);
```

This cfg can then be pass with a time series to cfc.estimate_cfc to perform the computation.
```matlab
cfc_results = cfc.estimate_cfc( signal,cfc_cfg );
```

Similarly, the whole comodulogram can be computed and plotted with:
```matlab
cmg_cfg = struct('sr',      sample_rate,...
    'lo_bounds',           [6 20],...
    'lo_bandwidth',        2,...
    'lo_step',             2,...
    'hi_bounds',           [20 150],...
    'hi_bandwidth',        'adaptive',...
    'hi_step',             10,...
    'metrics',             {{'MI'}});

cmg = cfc.estimate_comodulogram( signal,cmg_cfg );
cfc.plot.cmg(cmg)
```

There are examples of these and many more analyses in cfc.example using a simulated dataset, you can view and run through these examples by viewing the example script
```matlab
edit cfc.example
```
