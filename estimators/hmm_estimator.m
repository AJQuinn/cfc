function [ vpath ] = hmm_estimator(modulating_signal,modulated_signal)

signal = [modulating_signal, modulated_signal];
size(signal)

options = struct('K',2);
options.cyc = 100;
options.tol = 0.00001;
options.initcyc = 10;
options.order = 0;
options.timelag = 5;
options.embeddedlags = 100;
options.covtype = 'diag';

[hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist] = hmmmar (signal,size(signal,1),options);

vmat = zeros(size(vpath,1),options.K);
for i = 1:options.K
    vmat(:,i) = i - .3;
    vmat(vpath == i,i) = i + .3;
end
