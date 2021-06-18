
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

%% Run TrivialEstimate_2

%Set seed

rng(1996)

total_time = 25*40; %Therefore total time+1 total pieces of data, since I_0
...is at time = 0.

days = 0:total_time;

%Sinusoidal R_t

R_t = 1+0.1*cos(days*2*pi/25); 

w_s_o = [0 0.15 0.17 0.18 0.15 0.1 0.08 0.05 0.03 0.02 0.015 0.012...
    0.010 0.009 0.008 0.008 0.008]; %Original serial interval- gamma shaped
    
w_s_u = zeros(1, length(w_s_o));

w_s_u(end) = 1; %Dirac delta

w_s_t = zeros(total_time+1, length(w_s_o));

%%

w_s_o(1) = []; %Delete since our algorithm knows that probability of 
...serial=0 is 0.



tau = 8; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;

%Initial incidence
I_0 = 1000;
I = I_0;

%Generate incidence data

%Here we generate incidence data using the initial w_s_o only. Notice that
%we do not update the w_s value. This is because we want to predict how
%wrong our inference will be, namely that it will tend to having a lag of
%16 days. Meanwhile, we generate our serial interval as a weighted average
%(according to the time dissipated from when data was generated until the 
%end).

for t = 1:total_time
   
    I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_o)); 

%     I_new = R_t*Incidence_Generator_2(I,w_s); %Determinstic data
% %     generation
    I = [I, I_new];
    
    w_s_t(t+1, :) = (t*w_s_u +(total_time-t)*[0 w_s_o])/total_time;
    
end

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
        
        Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, w_s_t(k, :)); 
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

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1, ...
    'edgecolor', [1 1 1]);

hold on

h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');

%Plot Prior 95% CI. Useful when looking at failed epidemic

plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
h(3) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% plot([days(1) days(end)], [a*b a*b], 'b--')

%Plot ongoing 95% CI


%Labels

% title({['Trivial $R_t$ (=', num2str(R_t), ') inference with'];...
%     ['$w_s = [0$ ',num2str(w_s_o), '] \& $\tau=$ ', num2str(tau), ' days']})
% ylabel('$\bar{R}_t$')
% xlabel('Time, $t$ (days)')

%Legend

% legend(h([2]), '95 \% Confidence Interval')

legend(h([2 1 3]), '$\bar{R}_t$', '95 \% Posterior CI', '95 \% Prior CI', 'Location', 'Best')

Printer = 0;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Trivial_Estimate_CI.eps')

export_fig Trivial_Estimate_CI.eps -eps -r300 -painters -transparent

end

figure(2)
clf

%Plot Mean

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1, ...
    'edgecolor', [1 1 1]);

hold on

h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');

%Plot Prior 95% CI. Useful when looking at failed epidemic

plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
h(3) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% plot([days(1) days(end)], [a*b a*b], 'b--')

h(4) = plot(days, R_t, 'color', C(3, :))
h(4) = plot(tau:max(days)+tau, R_t, 'color', C(4, :))

%Plot ongoing 95% CI


%Labels

% title({['Trivial $R_t$ (=', num2str(R_t), ') inference with'];...
%     ['$w_s = [0$ ',num2str(w_s_o), '] \& $\tau=$ ', num2str(tau)', ' days']})
% ylabel('$\bar{R}_t$')
% xlabel('Time, $t$ (days)')

%Legend

% legend(h([2]), '95 \% Confidence Interval')

legend(h([2 1 3]), '$\bar{R}_t$', '95 \% Posterior CI', 'Location', 'Best')
axis([0 length(days) 0.25 1.75])

Printer = 0;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Trivial_Estimate_CI_Zoom.eps')

export_fig Trivial_Estimate_CI_Zoom.eps -eps -r300 -painters -transparent

end


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
saveas(gcf, 'Trivial_Estimate_Incidence')

export_fig Trivial_Estimate_Incidence.eps -eps -r300 -painters -transparent

end


%% Check that the mean serial time is line in growth- not that useful...
% Days_In_Serial = 0:length(w_s_o);
% 
% Weighted_Mean = w_s_t*Days_In_Serial';
% 
% figure(4)
% clf
% 
% plot(days, Weighted_Mean)
