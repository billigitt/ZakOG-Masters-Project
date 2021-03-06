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
Printer = 0;
%% 

[k, d, f] = Delta_Generator_Sigmoid_1(3, 0.065, 0.3);

w_o = [0.25 0.5 0.25];

N_o = length(w_o);

x = linspace(1,N_o,N_o);

w_a = d + w_o;

total_time = 30;

tau = 4;

w_control = w_o;

w_actual_1 = [w_o; w_a];

switch_behaviour = 16;

delay = 0;

update_behaviour = switch_behaviour + delay;

para_1_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_1_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_actual_1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial = struct('R_t', 0.8);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_correct, para_Trivial);

daysflip = [tau+1:total_time, total_time:-1:tau+1];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];

inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

%%
figure(1)

h(1) = plot(x, w_o, 'k');
hold on
h(2) = plot(x, w_a, 'color', C(1, :));

legend(h(1), 'Original')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

%%
figure(2)
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

legend(h([1 4 6 2 3]), '$R_t = 0.8$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)')

title('SI updates increasing mean by 0.065 and variance by 0.3')


%% R_t =2


[k, d, f] = Delta_Generator_Sigmoid_1(3, 0.065, 0.3);

w_o = [0.25 0.5 0.25];

N_o = length(w_o);

x = linspace(1,N_o,N_o);

w_a = d + w_o;

total_time = 30;

tau = 4;

w_control = w_o;

w_actual_1 = [w_o; w_a];

switch_behaviour = 16;

delay = 0;

update_behaviour = switch_behaviour + delay;

para_1_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_1_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_actual_1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial = struct('R_t', 2);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_correct, para_Trivial);

daysflip = [tau+1:total_time, total_time:-1:tau+1];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];

inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

%%
figure(3)

h(1) = plot(x, w_o, 'k');
hold on
h(2) = plot(x, w_a, 'color', C(1, :));

legend(h(1), 'Original')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

%%
figure(4)
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

title('SI updates increasing mean by 0.065 and variance by 0.3')


