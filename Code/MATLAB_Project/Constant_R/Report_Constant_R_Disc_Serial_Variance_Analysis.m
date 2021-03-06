%% Cleaning
clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 1.5)
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
C  = [0.3686 0.3098 0.6353; 0.2005 0.5593 0.7380; 0.4558 0.7897 0.6458;...
    0.8525 0.2654 0.3082; 0.6196 0.0039 0.2588];
addpath('../Functions', '~/Documents/MATLAB/export_fig', '~/Documents/MATLAB/matlab2tikz')
savepath ../Functions/pathdef.m
Printer = 1;
%% Report_Constant_R_Disc_Serial_Variance_Analysis

N_o = 15;

pd_o = makedist('Normal','mu',8,'sigma',2);

Length_Experiment = 1000;

Means = zeros(2, Length_Experiment);

for i = 1:Length_Experiment

    Means(:, i) = Report_Constant_R_Disc_Serial_Variance_Analyser(N_o, 8, 4, 1, pd_o, 100, 8, 40, 0, i);

end


%%
figure(1)
clf

h(1) = histogram(Means(1, :), 'BinEdges', [1.9:0.01:2.14], 'FaceColor', C(3, :), 'LineStyle', 'none');
xlabel('Final $R_t$ after SI change')
ylabel('Frequency')

hold on

h(2) = histogram(Means(2, :), 'BinEdges', [1.9:0.01:2.14], 'FaceColor', C(2, :), 'LineStyle', 'none');

legend(h([1 2]), '$\tilde{N}(8, 2)$ to $\tilde{N}(8, 4)$', '$\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

title({'Histograms of $R_t$ inference 60 days'; 'after SI change, without being updated'})

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_Variance_Analysis.eps -eps

end

function [Mean_After_Intervention] = Report_Constant_R_Disc_Serial_Variance_Analyser(N_o, mu, sigma_up, sigma_down, pd_o, Total_time, Tau, Switch_behaviour, Delay, Seed)

pd_a3 = makedist('Normal', 'mu', mu, 'sigma', sigma_up);
pd_a4 = makedist('Normal', 'mu', mu, 'sigma', sigma_down);

trunc_o = truncate(pd_o, 0, inf);

trunc_a3 = truncate(pd_a3, 0, inf);
trunc_a4 = truncate(pd_a4, 0, inf);

x = linspace(1,N_o,N_o);

w_s_o = pdf(trunc_o,x)/sum(pdf(trunc_o,x));
w_s_a3 = pdf(trunc_a3,x)/sum(pdf(trunc_a3,x));
w_s_a4 = pdf(trunc_a4,x)/sum(pdf(trunc_a4,x));

total_time = Total_time;

tau = Tau;

w_control = w_s_o;

w_actual_3 = [w_s_o; w_s_a3];

w_actual_4 = [w_s_o; w_s_a4];

switch_behaviour = Switch_behaviour;

delay = Delay;

update_behaviour = switch_behaviour + delay;

para_3_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_3, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_4_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_4, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[~, ~, ~, ~, ~, Mean_3_incorrect, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_3_incorrect, para_Trivial);

[~, ~, ~, ~, ~, Mean_4_incorrect, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_4_incorrect, para_Trivial);

Mean_3_Interest = mean(Mean_3_incorrect(end));
Mean_4_Interest = mean(Mean_4_incorrect(end));

Mean_After_Intervention = [Mean_3_Interest; Mean_4_Interest];

end
