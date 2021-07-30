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

%% Run Updating_SI_Constant_R_Truncated_SI

%We create a make-believe disease with a serial interval that suddenly
%truncates. We assume the initial serial interval is a Normal+ distribution
%with mean 5 and std 5.

N_o = 12;

pd = makedist('Normal','mu',5,'sigma',2);

trunc = truncate(pd, 0, inf);

x = linspace(1,N_o-1,N_o-1);

w_s_o = pdf(trunc,x)/sum(pdf(trunc,x));

% w_s_o = [0 w_s_o];

% writematrix(w_s_o,'test_serial.csv') %For comparison to EpiEstim

w_s_a = w_s_o;

w_s_a(floor(0.7*N_o):end) = 0;

w_s_a = w_s_a/sum(w_s_a);

tau = 8;

w_s_all_actual = [w_s_o; w_s_a];

w_s_all_recorded = [w_s_o; w_s_a];

switch_behaviour = 80;

delay = 0;

update_behaviour = switch_behaviour+ delay;

para = struct('seed', 209, 'total_time', 300, 'w_s_all_actual', w_s_all_actual, 'w_s_all_recorded', w_s_all_recorded, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[w_s_actual_Hybrid, w_s_recorded_Hybrid, ~, ~, ~, Mean_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Hybrid', para, para_Trivial);

[w_s_actual_Non_Hybrid, w_s_recorded_Non_Hybrid, I, Shape_Non_Hybrid, Scale_Non_Hybrid, Mean_Non_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para, para_Trivial);

writematrix([0 w_s_o],'Varying_Serial.csv')
writematrix(I','Varying_Serial_Data.csv')

%We note that (when generating data artifically) we want to always have 
%the data generated using hybrid SIs. Hybrid Generation is for this case
%but we note we only use it for one switch time.

[w_s_actual_Hybrid_Generation_only, w_s_recorded_Hybrid_Generation_only, ~, ~, ~, Mean_Hybrid_Generation_only, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Hybrid-Generation', para, para_Trivial);

figure(1)
clf

h(1) = plot(Mean_Non_Hybrid, 'k');
hold on
h(2) = plot(Mean_Hybrid, 'r');

h(3) = plot(Mean_Hybrid_Generation_only, 'b--');

% legend(h([1 2 3]), {'No Hybrid SIs', 'Hybrid SIs', 'Hybrid Generation'}, 'Location', 'NorthWest')

ylabel('$R_t$ inference')
xlabel('Time, $t$ (Days)')

axis([75 95 1.85 2.05])
I(2)
%%

para_Variable = struct('R_t', @(x) 2.0 - 0.5*(x>50));

[w_s_actual_Hybrid, w_s_recorded_Hybrid, ~, ~, ~, Mean_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Variable','Hybrid', para, para_Variable);

[w_s_actual_Non_Hybrid, w_s_recorded_Non_Hybrid, ~, ~, ~, Mean_Non_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para, para_Variable);

%We note that (when generating data artifically) we want to always have 
%the data generated using hybrid SIs. Hybrid Generation is for this case
%but we note we only use it for one switch time.

[w_s_actual_Hybrid_Generation_only, w_s_recorded_Hybrid_Generation_only, ~, ~, ~, Mean_Hybrid_Generation_only, ~, ~] = R_infer_disc_update_SI('Perfect', 'Variable','Hybrid-Generation', para, para_Variable);


figure(2)
clf

h(1) = plot(Mean_Non_Hybrid, 'k');
hold on
h(2) = plot(Mean_Hybrid, 'r');

h(3) = plot(Mean_Hybrid_Generation_only, 'b--');

legend(h([1 2 3]), {'No Hybrid SIs', 'Hybrid SIs', 'Hybrid Generation'}, 'Location', 'NorthEast')

ylabel('$R_t$ inference')
xlabel('Time, $t$ (Days)')
%%
% Go to https://shiny.dide.imperial.ac.uk/epiestim/ and upload data as well
% as downloading EstimatedR.csv to compare to shape and scale parameters
% through time. The curves should lie exactly on top of one another.

data = readtable('Varying_Serial_R.csv');

mu = data.Mean_R_;

sigma = data.Std_R_;

k = (mu.^2./sigma.^2);
theta = mu./k;
x = 0:0.001:6;
figure(5)
clf
plot(x, gampdf(x, k(150), theta(150)))

Printer = 0;


figure(6)
clf
h(1) = plot(k, 'g');
hold on
h(2) = plot(Shape_Non_Hybrid(tau+1:end), 'b--');
title('Shape parameter comparison')
ylabel('$k$')
xlabel('Time, $t$ (days)')
legend(h([1 2]), 'EpiEstim', 'My inference', 'Location', 'NorthWest')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Trivial_Estimate_k_Comparison')

export_fig Trivial_Estimate_k_Comparison.eps -eps -r300 -painters -transparent

end


figure(7)
clf
h(1) = plot(theta, 'g');
hold on
h(2) = plot(Scale_Non_Hybrid(tau+1:end), 'r--');
title('Scale parameter comparison')
ylabel('$\theta$')
xlabel('Time, $t$ (days)')
legend(h([1 2]), 'EpiEstim', 'My inference')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Trivial_Estimate_theta_Comparison')

export_fig Trivial_Estimate_theta_Comparison.eps -eps -r300 -painters -transparent

end

%Thoughts
% why is my one one longer than theirs?
% why do they use one less time pt for the summations (or appear to)?