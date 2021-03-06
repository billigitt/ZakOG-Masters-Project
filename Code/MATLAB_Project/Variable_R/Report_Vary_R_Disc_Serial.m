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
%% Report_Vary_R_Disc_Serial

%Here, we aim to create a large and comprehensive study looking at the
%effect of having good estimates on the SI, thinking about how the SI could
%change through time and thinking about the effect of linearly varying R
%values. We hope to find out how significant this change in R_t is with
%time.

%Plan

%Look at varying R without any change in SI.

%Vary R through the critical region R = 1

%Choose end points and then make it change through time depending on this

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

% Epidemiological and Inference parameters

total_time = 200;

tau = 8;

w_control = w_s_o;

w_actual_1 = [w_s_o; w_s_a1];

w_actual_2 = [w_s_o; w_s_a2];

w_actual_3 = [w_s_o; w_s_a3];

w_actual_4 = [w_s_o; w_s_a4];

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_control, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

Seed = 3;

para_1_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_2_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_3_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_3, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_4_incorrect = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_4, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_1_correct = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_actual_1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_2_correct = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_actual_2, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_3_correct = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_3, 'w_s_all_recorded', w_actual_3, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_4_correct = struct('seed', Seed, 'total_time', total_time, 'w_s_all_actual', w_actual_4, 'w_s_all_recorded', w_actual_4, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

R_start = 10;

R_end = 0.6;

days = 0:1:total_time;

para_Linear_Vary = struct('R_t',R_start + (R_end-R_start)*days/total_time);

[~, ~, ~, ~, ~, Mean_o, Upper_o, Lower_o] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_o, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_1_incorrect, para_Linear_Vary);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_1_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_2_incorrect, Upper_2_incorrect, Lower_2_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_2_incorrect, para_Linear_Vary);

[~, ~, I_2, ~, ~, Mean_2_correct, Upper_2_correct, Lower_2_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_2_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_3_incorrect, Upper_3_incorrect, Lower_3_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_3_incorrect, para_Linear_Vary);

[~, ~, I_3, ~, ~, Mean_3_correct, Upper_3_correct, Lower_3_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_3_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_4_incorrect, Upper_4_incorrect, Lower_4_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_4_incorrect, para_Linear_Vary);

[~, ~, I_4, ~, ~, Mean_4_correct, Upper_4_correct, Lower_4_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_4_correct, para_Linear_Vary);


daysflip = [tau+1:total_time, total_time:-1:tau+1];

inBetween_o = [Lower_o(tau+1:total_time), fliplr(Upper_o(tau+1:total_time))];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];
inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

inBetween_2_incorrect = [Lower_2_incorrect(tau+1:total_time), fliplr(Upper_2_incorrect(tau+1:total_time))];
inBetween_2_correct = [Lower_2_correct(tau+1:total_time), fliplr(Upper_2_correct(tau+1:total_time))];

inBetween_3_incorrect = [Lower_3_incorrect(tau+1:total_time), fliplr(Upper_3_incorrect(tau+1:total_time))];
inBetween_3_correct = [Lower_3_correct(tau+1:total_time), fliplr(Upper_3_correct(tau+1:total_time))];

inBetween_4_incorrect = [Lower_4_incorrect(tau+1:total_time), fliplr(Upper_4_incorrect(tau+1:total_time))];
inBetween_4_correct = [Lower_4_correct(tau+1:total_time), fliplr(Upper_4_correct(tau+1:total_time))];


%

figure(1)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_o, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_o(tau+1:total_time), 'k');
h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 2 4]), 'True $R_t$', 'Mean $\tilde{R}_t$', '95 \% CI', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_o.eps -eps


end


%
figure(2)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);hold on
hold on
h(2) = fill(daysflip, inBetween_1_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_1_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_1_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_1_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_1_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(10, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_1.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

xlabel('')
ylabel('')

figure(3)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);hold on
hold on
h(2) = fill(daysflip, inBetween_2_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_2_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_2_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_2_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_2_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(6, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_2.eps -eps -opengl

end


figure(4)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);hold on
hold on
h(2) = fill(daysflip, inBetween_3_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_3_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_3_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_3_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_3_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 4)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_3.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

xlabel('')
ylabel('')

figure(5)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);hold on

hold on

h(2) = fill(daysflip, inBetween_4_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_4_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_4_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_4_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_4_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_4.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end
ylabel('')
%

figure(6)

h(1) = plot(x, w_s_o, 'k');
hold on
h(2) = plot(x, w_s_a1, 'color', C(1, :));
h(3) = plot(x, w_s_a2, 'color', C(1, :), 'LineStyle', '--');
h(4) = plot(x, w_s_a3, 'color', C(1, :), 'LineStyle', '-.');
h(5) = plot(x, w_s_a4, 'color', C(1, :), 'LineStyle', ':');

axis([0 15 0 0.45])

legend(h(1), 'Original', 'Location', 'best')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_0.eps -eps

% 
% matlab2tikz('triagle.tex'); % black background

end



% Now try varying R the opposite way
%
%
%
%
%
%
%
%
%
%
%


R_start = 0.6;

R_end = 10;

days = 0:1:total_time;

para_Linear_Vary = struct('R_t',R_start + (R_end-R_start)*days/total_time);

[~, ~, ~, ~, ~, Mean_o, Upper_o, Lower_o] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_o, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_1_incorrect, para_Linear_Vary);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_1_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_2_incorrect, Upper_2_incorrect, Lower_2_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_2_incorrect, para_Linear_Vary);

[~, ~, I_2, ~, ~, Mean_2_correct, Upper_2_correct, Lower_2_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_2_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_3_incorrect, Upper_3_incorrect, Lower_3_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_3_incorrect, para_Linear_Vary);

[~, ~, I_3, ~, ~, Mean_3_correct, Upper_3_correct, Lower_3_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_3_correct, para_Linear_Vary);

[~, ~, ~, ~, ~, Mean_4_incorrect, Upper_4_incorrect, Lower_4_incorrect] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_4_incorrect, para_Linear_Vary);

[~, ~, I_4, ~, ~, Mean_4_correct, Upper_4_correct, Lower_4_correct] = R_infer_disc_update_SI('Perfect', 'Variable','Non-Hybrid', para_4_correct, para_Linear_Vary);


daysflip = [tau+1:total_time, total_time:-1:tau+1];

inBetween_o = [Lower_o(tau+1:total_time), fliplr(Upper_o(tau+1:total_time))];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];
inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

inBetween_2_incorrect = [Lower_2_incorrect(tau+1:total_time), fliplr(Upper_2_incorrect(tau+1:total_time))];
inBetween_2_correct = [Lower_2_correct(tau+1:total_time), fliplr(Upper_2_correct(tau+1:total_time))];

inBetween_3_incorrect = [Lower_3_incorrect(tau+1:total_time), fliplr(Upper_3_incorrect(tau+1:total_time))];
inBetween_3_correct = [Lower_3_correct(tau+1:total_time), fliplr(Upper_3_correct(tau+1:total_time))];

inBetween_4_incorrect = [Lower_4_incorrect(tau+1:total_time), fliplr(Upper_4_incorrect(tau+1:total_time))];
inBetween_4_correct = [Lower_4_correct(tau+1:total_time), fliplr(Upper_4_correct(tau+1:total_time))];



figure(7)

h(1) = plot(x, w_s_o, 'k');
hold on
h(2) = plot(x, w_s_a1, 'color', C(1, :));
h(3) = plot(x, w_s_a2, 'color', C(1, :), 'LineStyle', '--');
h(4) = plot(x, w_s_a3, 'color', C(1, :), 'LineStyle', '-.');
h(5) = plot(x, w_s_a4, 'color', C(1, :), 'LineStyle', ':');

axis([0 15 0 0.45])

legend(h(1), 'Original', 'Location', 'best')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_0.eps -eps

% 
% matlab2tikz('triagle.tex'); % black background

end

figure(8)
clf
h(8) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_1_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_1_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_1_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_1_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_1_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([8 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(10, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_1.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

xlabel('')
ylabel('')

figure(9)
clf
h(8) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_2_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_2_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_2_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_2_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_2_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([8 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(6, 2)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_2.eps -eps -opengl

end

figure(10)
clf
h(8) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_3_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = fill(daysflip, inBetween_3_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_3_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_3_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_3_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([8 4 6 2 3]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 4)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_3.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

xlabel('')
ylabel('')

figure(11)
clf
% h(1) = plot([0 total_time], [para_Linear_Vary.R_t para_Linear_Vary.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on

h(1) = fill(daysflip, inBetween_4_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(2) = fill(daysflip, inBetween_4_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(3) = plot(tau+1:total_time, Mean_4_incorrect(tau+1:total_time), 'color', C(4, :));
h(4) = plot(tau+1:switch_behaviour, Mean_4_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(5) = plot(switch_behaviour:total_time, Mean_4_correct(switch_behaviour:total_time), 'color', C(1, :));

h(6) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

h(7) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);


xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([7 3 5 1 2]), 'True $R_t$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', '95 \% CI (original SI)', '95 \% CI (updated SI)', 'Location', 'best')

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_4.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

ylabel('')

figure(12)
clf
h(1) = plot(0:1:total_time, para_Linear_Vary.R_t, 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_o, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_o(tau+1:total_time), 'k');
h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$')

legend(h([1 4 2]), 'True $R_t$', 'Mean $\tilde{R}_t$', '95 \% CI', 'Location', 'best');

title('SI updates from $\tilde{N}(8, 2)$ to $\tilde{N}(8, 1)$')

if Printer == 1

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_o.eps -eps


end

%
perm  = [2 5 3 6];

figure(13)
clf
ax = zeros(6, 1);
for i = 1:6
    ax(i) = subplot(2, 3, i);
    
    if i == 4
       
        axis off
       
    end
end
% Now copy contents of each figure over to destination figure
% Modify position of each axes as it is transferred
figure(1)
    g = get(gcf,'Children');
    g(1) = []; %Removes legend
    newh = copyobj(g,13);
    for j = 1:length(newh)
posnewh = get(newh(j),'Position');
possub  = get(ax(1),'Position');
set(newh(j),'Position',...
[possub(1) possub(2) possub(3) possub(4)])


    end
    delete(ax(1));

for i = 2:5
    figure(i)
    g = get(gcf,'Children');
    
    if i ~= 5
        g(1) = []; %Removes legend
        newh = copyobj(g,13);
        
        posnewh = get(newh(1),'Position');
    possub  = get(ax(perm(i-1)),'Position');
    set(newh(1),'Position',...
    [possub(1) possub(2) possub(3) possub(4)])
        
    else
        
        newh = copyobj(g,13);
        
        set(newh(1),'Position',...
     [0.13 0.267 0.175 0.175])
        
        posnewh = get(newh(2),'Position');
    possub  = get(ax(perm(i-1)),'Position');
    set(newh(2),'Position',...
    [possub(1) possub(2) possub(3) possub(4)])
        
    end
    
    delete(ax(perm(i-1)));
end
figure(13)
set(gcf, 'Units', 'centimeters', 'Position', [0 0 37 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);

if Printer == 5

% Save figure
export_fig Report_Vary_Lin_Decrease_R_Disc_Serial_All.eps -eps -opengl

end


%

perm  = [2 5 3 6];

figure(14)
clf
ax = zeros(6, 1);
for i = 1:6
    ax(i) = subplot(2, 3, i);
    
    if i == 4
       
        axis off
       
    end
end
% Now copy contents of each figure over to destination figure
% Modify position of each axes as it is transferred
figure(12)
    g = get(gcf,'Children');
    g(1) = []; %Removes legend
    newh = copyobj(g,14);
    for j = 1:length(newh)
posnewh = get(newh(j),'Position');
possub  = get(ax(1),'Position');
set(newh(j),'Position',...
[possub(1) possub(2) possub(3) possub(4)])


    end
    delete(ax(1));

for i = 8:11
    figure(i)
    g = get(gcf,'Children');
    
    if i ~= 11
        g(1) = []; %Removes legend
        newh = copyobj(g,14);
        
        posnewh = get(newh(1),'Position');
    possub  = get(ax(perm(i-7)),'Position');
    set(newh(1),'Position',...
    [possub(1) possub(2) possub(3) possub(4)])
        
    else
        
        newh = copyobj(g,14);
        
        set(newh(1),'Position',...
     [0.13 0.267 0.175 0.175])
        
        posnewh = get(newh(2),'Position');
    possub  = get(ax(perm(i-7)),'Position');
    set(newh(2),'Position',...
    [possub(1) possub(2) possub(3) possub(4)])
        
    end
    
    delete(ax(perm(i-7)));
end
figure(14)
set(gcf, 'Units', 'centimeters', 'Position', [0 0 37 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);

if Printer == 0

% Save figure
export_fig Report_Vary_Lin_Increase_R_Disc_Serial_All.eps -eps -opengl

end