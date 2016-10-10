% speedup
clear all; clc; close all;

% #threads
processors = [1 2 4 8];

% Experimental data

% time; 1 processor; 2 processor; 4 processor; processor 8: processor 16;
time = [2.4455799999996088 1.4056020000134595 0.8547910000197589 0.5539709999575280;
        2.4428740000003017 1.4161969999549910 0.8663079999969341 0.5506129999994300;
        2.4480209999601357 1.4047319999663159 0.8514680000371300 0.5663649999769405;
        2.4452739999978803 1.4123290000134148 0.8563059999723919 0.6018219999969006;
        2.4454329999862239 1.4169799999799579 0.8768750000162981 0.6547780000255443];  
time_mean = mean(time);

% Calculating speedup, S = Tinit/Tnew
speedup = time_mean(1)./time_mean(:);

figure
plot(processors,time_mean,processors,speedup, 'LineWidth',2)
set(gca,'FontSize',13, 'FontWeight', 'bold');
title(['Time cost and speedup for Quick sort in Parallel']);
ylabel('Speedup and time cost') % label left y-axis
xlabel('Processors') % label x-axis
legend('Time cost','Speedup')
grid on
