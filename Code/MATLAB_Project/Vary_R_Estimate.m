
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
%% Run Vary_R_Estimate

%Set seed

rng(1996)

total_time = 80; %Therefore total time+1 total pieces of data, since I_0
...is at time = 0.

days = 0:total_time;


%Variable R_t

R_t = 0.75 + 3*exp(-0.25*days);

w_s = [0 0.1 0.2 0.3 0.2 0.1 0.05 0.03 0.02]; %Serial interval, e.g. odds
...of infecting after 1 days is 1/3.

writematrix(w_s,'Trivial_Serial.csv') %For comparison to EpiEstim

w_s(1) = []; %Delete since our algorithm knows that probability of 
...serial=0 is 0.



tau = 8; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;

%Initial incidence
I_0 = 10;
I = I_0;

%Generate incidence data

%Incidence_Generator_2 takes the whole time series data for I, extracts the
%most recent (last 8 days in this case) incidence reports and then computes
%the dot product with the serial. [The serial is reversed before taking the
%dot product]

for t = 1:total_time
   
    I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s));

%     I_new = R_t*Incidence_Generator_2(I,w_s); %Determinstic data
% %     generation
    I = [I, I_new];
    
end

writematrix(I','Trivial_Data.csv') %For comparison to EpiEstim

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
        
        Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s]); 
        %Takes I_k, I_{k-1} down to length of Serial
        %The 0 is included again. Previously it was ommitted, since then we
        %were *generating* the *following* day so we wanted the 'current'
        %(but really previous day) to be included. Here we are actually
        %thinking about the current day, i.e. k and so we re-include the 0,
        %to illustrate that the serial must be 0 at this value.
    end
    
    I(1) = I_0; %First case is imported
    
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

h(1) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');
hold on



daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
h(2) = fill(daysflip, inBetween, 'k', 'FaceAlpha', 0.25, ...
    'edgealpha', 0);
h(3) = plot(days, R_t, 'color', C(4, :));
%Labels

plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% plot([days(1) days(end)], [a*b a*b], 'b--')

title({['$R_t$ inference with'];...
    ['$w_s = [$',num2str(w_s), ']']})
ylabel('$\bar{R}_t$')
xlabel('Time, $t$ (days)')

%Legend

legend(h([1 2 3]), '$\bar{R}_t$', '95\% confidence interval', 'True $R_t(t) = 0.75 + 3e^{-0.25t}$', 'Location', 'West')

Printer = 1;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Variable_Estimate_CI.pdf')

export_fig Variable_Estimate_CI.pdf -pdf -r1500 -opengl

end

figure(2)
clf
%Plot Mean

h(1) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');
hold on

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
h(2) = fill(daysflip, inBetween, 'k', 'FaceAlpha', 0.25, ...
    'edgealpha', 0);
h(3) = plot(days, R_t, 'color', C(4, :));
%Labels

plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% plot([days(1) days(end)], [a*b a*b], 'b--')

title({['$R_t$ inference with'];...
    ['$w_s = [$',num2str(w_s), ']']})
ylabel('$\bar{R}_t$')
xlabel('Time, $t$ (days)')

%Legend

legend(h([1 2 3]), '$\bar{R}_t$', '95\% confidence interval', 'True $R_t(t) = 0.75 + 3e^{-0.25t}$', 'Location', 'North')

axis([0 60 0 4])

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Variable_Estimate_CI_Zoom.pdf')

export_fig Variable_Estimate_CI_Zoom.pdf -pdf -r1500 -opengl

end

%%

figure(3)
clf
h(1) = bar(days, I, 'FaceColor', [.5 .5 .5]);
% hold on
% plot(days, R_t.^days, 'r')
title('Simulated incidence data')
ylabel('Incidence')
xlabel('Time, $t$ (days)')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Variable_Estimate_Incidence.eps')

export_fig Variable_Estimate_Incidence.eps -eps -r300 -painters -transparent

end

figure(4)
plot(0:8, [0 w_s], 'k')
title('Serial Interval (zero uncertainty)')
ylabel('Probability')
xlabel('Time, $t$ (days)')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Variable_Estimate_Certain_Serial.eps')

export_fig Variable_Estimate_Certain_Serial.eps -eps -r300 -painters -transparent

end

%% Comparison to EpiEstim

% Go to https://shiny.dide.imperial.ac.uk/epiestim/ and upload data as well
% as downloading EstimatedR.csv to compare to shape and scale parameters
% through time. The curves should lie exactly on top of one another.

data = readtable('Estimated_R_vary.csv');

mu = data.Mean_R_;

sigma = data.Std_R_;

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
saveas(gcf, 'Variable_Estimate_parameter_Comparison.eps')

export_fig Variable_Estimate_parameter_Comparison.eps -eps -r300 -painters -transparent

end