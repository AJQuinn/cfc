function plot_spec(signal,low,actual)

figure('Position',[1,1,500,800]);
axes('position',[.1,.35,.8,.6])

spectrogram(signal,56);
view(-90, 90)

axes('position',[.1,.1,.8,.2])
plot(low);
hold on;
plot(actual,'r');
axis tight;set(gca,'Xdir','reverse');
