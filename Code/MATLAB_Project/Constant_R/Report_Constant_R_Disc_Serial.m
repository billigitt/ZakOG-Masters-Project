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
Printer = 2;
%% Report_Constant_R_Disc_Serial

%Here, we aim to create a large and comprehensive study looking at the
%effect of having good estimates on the SI, thinking about how the SI could
%change through time and thinking about the effect of different constant R
%values. We hope to find out how significant this change in R_t is with
%time.

%Plan

%R = 0.8, 1, 1.1, 2, (5?)
% Change SI by truncation but also by just reducing the mean and variance.
% Perhaps twice.

%% 2 SIs, shift SI up and down, compare good and bad SIs, R= 0.8, 1, 1.1, 2

%% Create SIs

N_o = 15;

pd_o = makedist('Normal','mu',8,'sigma',2);
pd_a1 = makedist('Normal','mu',10,'sigma',2);
pd_a2 = makedist('Normal','mu',6,'sigma',2);
pd_a3 = makedist('Normal', 'mu', 8, 'sigma', 4);
pd_a4 = makedist('Normal', 'mu', 8, 'sigma', 1);

trunc_o = truncate(pd_o, 0, inf);
trunc_a1 = truncate(pd_a1, 0, inf);
trunc_a2 = truncate(pd_a2, 0, inf);
trunc_a3 = truncate(pd_a3, 0, inf);
trunc_a4 = truncate(pd_a4, 0, inf);

x = linspace(1,N_o,N_o);

w_s_o = pdf(trunc_o,x)/sum(pdf(trunc_o,x));
w_s_a1 = pdf(trunc_a1,x)/sum(pdf(trunc_a1,x));
w_s_a2 = pdf(trunc_a2,x)/sum(pdf(trunc_a2,x));
w_s_a3 = pdf(trunc_a3,x)/sum(pdf(trunc_a3,x));
w_s_a4 = pdf(trunc_a4,x)/sum(pdf(trunc_a4,x));

%% Epidemiological and Inference parameters

total_time = 100;

tau = 8;

w_control = w_s_o;

w_actual_1 = [w_s_o; w_s_a1];

w_actual_2 = [w_s_o; w_s_a2];

w_actual_3 = [w_s_o; w_s_a3];

w_actual_4 = [w_s_o; w_s_a4];

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

para_1_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_2_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_3_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_3, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_4_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_4, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_1_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_actual_1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_2_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_actual_2, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_3_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_3, 'w_s_all_recorded', w_actual_3, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_4_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_4, 'w_s_all_recorded', w_actual_4, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[~, ~, ~, ~, ~, Mean_0_default, Upper_0_default, Lower_0_default] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_correct, para_Trivial);

[~, ~, ~, ~, ~, Mean_2_incorrect, Upper_2_incorrect, Lower_2_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_2_incorrect, para_Trivial);

[~, ~, I_2, ~, ~, Mean_2_correct, Upper_2_correct, Lower_2_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_2_correct, para_Trivial);

[~, ~, ~, ~, ~, Mean_3_incorrect, Upper_3_incorrect, Lower_3_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_3_incorrect, para_Trivial);

[~, ~, I_3, ~, ~, Mean_3_correct, Upper_3_correct, Lower_3_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_3_correct, para_Trivial);

[~, ~, ~, ~, ~, Mean_4_incorrect, Upper_4_incorrect, Lower_4_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_4_incorrect, para_Trivial);

[~, ~, I_4, ~, ~, Mean_4_correct, Upper_4_correct, Lower_4_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_4_correct, para_Trivial);


daysflip = [tau+1:total_time, total_time:-1:tau+1];

inBetween_0_default = [Lower_0_default(tau+1:total_time), fliplr(Upper_0_default(tau+1:total_time))];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];
inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

inBetween_2_incorrect = [Lower_2_incorrect(tau+1:total_time), fliplr(Upper_2_incorrect(tau+1:total_time))];
inBetween_2_correct = [Lower_2_correct(tau+1:total_time), fliplr(Upper_2_correct(tau+1:total_time))];

inBetween_3_incorrect = [Lower_3_incorrect(tau+1:total_time), fliplr(Upper_3_incorrect(tau+1:total_time))];
inBetween_3_correct = [Lower_3_correct(tau+1:total_time), fliplr(Upper_3_correct(tau+1:total_time))];

inBetween_4_incorrect = [Lower_4_incorrect(tau+1:total_time), fliplr(Upper_4_incorrect(tau+1:total_time))];
inBetween_4_correct = [Lower_4_correct(tau+1:total_time), fliplr(Upper_4_correct(tau+1:total_time))];


%%
figure(1)

h(1) = plot(x, w_s_o, 'k');
hold on
h(2) = plot(x, w_s_a1, 'color', C(1, :));
h(3) = plot(x, w_s_a2, 'color', C(1, :), 'LineStyle', '--');
h(4) = plot(x, w_s_a3, 'color', C(1, :), 'LineStyle', '-.');
h(5) = plot(x, w_s_a4, 'color', C(1, :), 'LineStyle', ':');

axis([0 15 0 0.45])

legend(h(1), 'Original')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_0.eps -eps

% 
% matlab2tikz('triagle.tex'); % black background

end
%%
figure(2)
clf
h(1) = plot([0 total_time], [para_Trivial.R_t para_Trivial.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_0_default, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_0_default(tau+1:total_time), 'color', 'k');

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 2]), '$R_t = 2$', 'Mean $\tilde{R}_t$', '95 \% CI')

title('Standard $R_t$ Inference')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_0_Default.eps -eps -opengl

end

%%
figure(3)
clf
h(1) = plot([0 total_time], [para_Trivial.R_t para_Trivial.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_1_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_1_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_1_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_1_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_1_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), '$R_t = 2$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(10, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_1.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

figure(4)
clf
h(1) = plot([0 total_time], [para_Trivial.R_t para_Trivial.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_2_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_2_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_2_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_2_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_2_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), '$R_t = 2$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(6, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_2.eps -eps -opengl

end

figure(5)
clf
h(1) = plot([0 total_time], [para_Trivial.R_t para_Trivial.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_3_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_3_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_3_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_3_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_3_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), '$R_t = 2$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 4)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_3.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

figure(6)
clf
h(1) = plot([0 total_time], [para_Trivial.R_t para_Trivial.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_4_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_4_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_4_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_4_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_4_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), '$R_t = 2$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_4.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

%%

