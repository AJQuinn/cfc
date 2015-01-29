function plot_pacman_sw(pac_obj)

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
subplot(422)
plot(pac_obj.results.time_vect,pac_obj.results.plv);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('PLV');
subplot(423)
plot(pac_obj.results.time_vect,pac_obj.results.esc.^2);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('ESC');
subplot(424);
plot(pac_obj.results.time_vect,pac_obj.results.nesc.^2);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('NESC');
subplot(425)
plot(pac_obj.results.time_vect,pac_obj.results.mi);hold on;
mmi = max(pac_obj.results.mi);
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse.*(1.2*mmi) ,'r')
title('Modulation Index')
subplot(426)
plot(pac_obj.results.time_vect,pac_obj.results.voytek);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
title('AEC')
subplot(427)
plot(pac_obj.results.time_vect,pac_obj.results.aec);hold on;
plot(pac_obj.results.time_vect,pac_obj.results.true_timecourse,'r')
