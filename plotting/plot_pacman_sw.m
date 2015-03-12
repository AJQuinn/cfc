function plot_pacman_sw(pac_obj, perms)

figure;
subplot(321)
plot(pac_obj.signals.time_vect,pac_obj.signals.theta);
title('Theta');
subplot(322)
plot(pac_obj.signals.time_vect,pac_obj.signals.gamma)
title('Gamma');
subplot(323)
plot(pac_obj.signals.time_vect,pac_obj.signals.theta_phase);
title('Theta Phase');
subplot(324)
plot(pac_obj.signals.time_vect,pac_obj.signals.gamma_amp);
title('Gamma Amp');
subplot(325)
plot(pac_obj.signals.time_vect,pac_obj.signals.gamma_amp_phase);
title('Gamma Amp Phase');
subplot(326)
plot(pac_obj.signals.time_vect,pac_obj.signals.gamma_amp_theta);
title('Gamma Amp Theta');

figure;
subplot(421)
plot(pac_obj.results.time_vect,pac_obj.results.glm');hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('GLM');
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.glm_cv, perms.glm_cv]);
end

subplot(422)
plot(pac_obj.results.time_vect,pac_obj.results.plv);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('PLV');
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.plv_cv, perms.plv_cv]);
end

subplot(423)
plot(pac_obj.results.time_vect,pac_obj.results.esc);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('Envelope Signal Corr');
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.esc_cv, perms.esc_cv]);
end

subplot(424);
plot(pac_obj.results.time_vect,pac_obj.results.nesc);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('Normalised Envelope Signal Corr');
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.nesc_cv, perms.nesc_cv]);
end

subplot(425)
plot(pac_obj.results.time_vect,pac_obj.results.mi);hold on;
mmi = max(pac_obj.results.mi);
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse.*(1.2*mmi) ,'r')
title('Modulation Index')
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.mi_cv, perms.mi_cv]);
end

subplot(426)
plot(pac_obj.results.time_vect,pac_obj.results.voytek);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('Voytek')
subplot(427)
plot(pac_obj.results.time_vect,pac_obj.results.aec);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('Amplitude Envelope Corr')
if nargin == 2
    line([pac_obj.results.time_vect(1),pac_obj.results.time_vect(end)],...
         [perms.aec_cv, perms.aec_cv]);
end



