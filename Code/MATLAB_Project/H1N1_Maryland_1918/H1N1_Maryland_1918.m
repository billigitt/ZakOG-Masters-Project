%% Cleaning
clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 2)
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
C  = [0.3686 0.3098 0.6353; 0.2005 0.5593 0.7380; 0.4558 0.7897 0.6458;...
    0.8525 0.2654 0.3082; 0.6196 0.0039 0.2588];
addpath('../Functions', '~/Documents/MATLAB/export_fig')  
savepath ../Functions/pathdef.m

%% Run H1N1_Maryland_1918

%Read data

inference_data = readtable('Incidence_H1N1_Maryland_1918.csv');
serial_data = readtable('Serial_Interval_H1N1_Maryland_1918.csv');

local = inference_data.local;
imported = inference_data.imported;
total_time = length(inference_data.local)-1; % -1 because we define the start from

I = imported'+local'; %Prelimary calculation. NEED TO CHANGE.
%time 0

days = 0:total_time;

w_s = (serial_data.x)'; %Serial interval, e.g. odds
...of infecting after 1 days is 1/3.

tau = 7; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;



%Generate incidence data

Shape = zeros(1, total_time+1);
Scale = zeros(1, total_time+1);
Mean = zeros(1, total_time+1);
Upper = zeros(1, total_time+1);
Lower = zeros(1, total_time+1);

for t = tau+1:total_time+1 %tau+1 (Day tau) is the first you can start 
    
    ...from. t is the index, NOT the day.
%     disp(I(t-tau:t))
    Shape(t) = a + sum(I(t-tau+1:t)); % IS +1 RIGHT?? shape param at time t = i
    Scale(t) = 0;
    
    %Calculate summation of Lambdas
    
    for k  = t-tau+1:t %Support material suggests that we only look at 
        ...these values, contrary to Thompson-Stockwin paper
        
        I_relevant = I(1:k); %definitely want from k down, may want all
        ...the way down to 1
        
        Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, w_s); 
        %Takes I_k, I_{k-1} down to length of Serial
        %The 0 is included again. Previously it was ommitted, since then we
        %were *generating* the *following* day so we wanted the 'current'
        %(but really previous day) to be included. Here we are actually
        %thinking about the current day, i.e. k and so we re-include the 0,
        %to illustrate that the serial must be 0 at this value.
    end
    
%     I(1) = I_0; %First case is imported
    
    Scale(t) = 1/(Scale(t)+(1/b));
    
    Mean(t) = Scale(t)*Shape(t); %Mean of gamma distn is product of two 
    ...pars
    
    Upper(t) = gaminv(0.975, Shape(t), Scale(t));
    
    Lower(t) = gaminv(0.025, Shape(t), Scale(t));
    
end

%% Plots

figure(1)
clf

%Plot Mean

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1,'edgecolor', [1 1 1]);

hold on

h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');

h(3) = plot([tau total_time], [1 1], 'k--', 'LineWidth', 1);

%Labels

% plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
% h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% % plot([days(1) days(end)], [a*b a*b], 'b--')

title({'Real-time $R_t$ inference of 1918 H1N1';['outbreak in Maryland, USA ($\tau =$ ', num2str(tau), ' days$)$']})
ylabel('$\bar{R}_t$')
xlabel('Time, $t$ (days)')

%Legend

legend(h([1 2]), '$\bar{R}_t$', '95\% Posterior CI', 'Location', 'Best')

Printer = 0;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'H1N1_Maryland_1918_Estimate_CI.eps')

export_fig H1N1_Maryland_1918_Estimate_CI.eps -eps -r300 -painters -transparent

end

figure(2)
clf

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];

yyaxis left

h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1,'edgecolor', [1 1 1]);

hold on

h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k', 'LineStyle', '-');

h(3) = plot([tau total_time], [1 1], 'k--', 'LineWidth', 1);

title({'Real-time $R_t$ inference of 1918 H1N1';['outbreak in Maryland, USA ($\tau =$ ', num2str(tau), ' days$)$']})
ylabel('$\tilde{R}_t$')
xlabel('Time, $t$ (days)')

yyaxis right

h(4) = bar(days, I, 'FaceColor', C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.3);

ylabel('Incidence')

ax = gca;
ax.YAxis(1).Color = 'k';

ax.YAxis(2).Color = C(1, :);

%Labels

% plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
% h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% % plot([days(1) days(end)], [a*b a*b], 'b--')



%Legend

legend(h([1 2 4]), '$\tilde{R}_t$', '95\% Posterior CI', 'Incidence Data', 'Location', 'Best')
% axis([0 total_time+tau 0 max(Upper)])

if Printer == 0
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'H1N1_Maryland_1918_Estimate_CI_Zoom.eps')

export_fig H1N1_Maryland_1918_Estimate_CI_Zoom.eps -eps -r300 -opengl

end

%%

figure(3)
clf
h(1) = bar(days, I, 'FaceColor', [.5 .5 .5]);
% hold on
% plot(days, R_t.^days, 'r')
title('Incidence data for H1N1 outbreak in Maryland, USA')
ylabel('Incidence')
xlabel('Time, $t$ (days)')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'H1N1_Maryland_1918_Estimate_Incidence.eps')

export_fig H1N1_Maryland_1918_Estimate_Incidence.eps -eps -r300 -painters -transparent

end

ix = find(w_s~=0, 1, 'last');

W_s = w_s(1:ix);

figure(4)
plot(0:length(W_s)-1, W_s, 'k')
title('EpiEstim Serial Interval')
ylabel('Probability')
xlabel('Time, $t$ (days)')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'H1N1_Maryland_1918_Estimate_Certain_Serial.eps')

export_fig H1N1_Maryland_1918_Estimate_Certain_Serial.eps -eps -r300 -painters -transparent

end

%% Comparison to EpiEstim

% Go to https://shiny.dide.imperial.ac.uk/epiestim/ and upload data as well
% as downloading EstimatedR.csv to compare to shape and scale parameters
% through time. The curves should lie exactly on top of one another.

inference_data = readtable('Estimated_R_H1N1_Maryland_1918.csv');

mu = inference_data.Mean_R_;

sigma = inference_data.Std_R_;

k = (mu.^2./sigma.^2);
theta = mu./k;
% x = 0:0.001:6;
% figure(5)
% clf
% plot(x, gampdf(x, k(end), theta(end)))

Printer = 1;


figure(5)
clf
subplot(1, 2, 1)
h(1) = plot(k,'color', C(3, :));
hold on
h(2) = plot(Shape(tau+1:end), 'k', 'LineStyle', '--');
title('Shape parameter comparison')
ylabel('$k$')
xlabel('Time, $t$ (days)')
legend(h([1 2]), 'EpiEstim', 'My inference', 'Location', 'NorthEast')

subplot(1, 2, 2)
h(1) = plot(theta, 'color', C(3, :));
hold on
h(2) = plot(Scale(tau+1:end), 'k', 'LineStyle', '--');
title('Scale parameter comparison')
ylabel('$\theta$')
xlabel('Time, $t$ (days)')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 30 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'H1N1_Maryland_1918_Estimate_parameter_Comparison.eps')

export_fig H1N1_Maryland_1918_Estimate_parameter_Comparison.eps -eps -r300 -painters -transparent

end