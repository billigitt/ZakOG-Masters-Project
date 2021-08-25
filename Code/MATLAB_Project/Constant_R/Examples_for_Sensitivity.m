%% Cleaning
clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 3)
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

N_o = 10;


N=10;

M=2*N+1;

B = [N+1 2; 2*N+1 3];

A = blkdiag(B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B);

mu = linspace(1, -1, 2*N+1);

x = zeros(4*N+2, 1);

for i = 1:2*N+1
   
    x(2*i-1) = 0;
    
    x(2*i) = 6*mu(i)/(N*(N+1));
    
end

phi = A\x;

a = phi(1:2:4*N+1);

b = phi(2:2:4*N+2);

Delta_Mat = repmat(b, 1, N) + a.*(1:N);

Delta_Mat = reshape(Delta_Mat', [1 N 2*N+1]);

w_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

w_a_Matrix = repmat(w_o, 1, 1, 2*N+1)+ Delta_Mat;

w_all_actual = zeros(2, N, 2*N+1);

w_all_actual(1, :, :) = repmat(w_o, 1, 1, 2*N+1);

w_all_actual(2, :, :) = w_a_Matrix;

total_time = 100;

tau = 7;

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

I_0 = 100000;

para = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_all_actual, 'w_s_all_recorded', w_all_actual, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);

w_s_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

w_s_a1 = w_o + Delta_Mat(1, :, 4);


w_s_a2 = w_o + Delta_Mat(1, :, 18);

x = linspace(1,N_o,N_o);
%% Epidemiological and Inference parameters

total_time = 100;

tau = 8;

w_control = w_s_o;

w_actual_1 = [w_s_o; w_s_a1];

w_actual_2 = [w_s_o; w_s_a2];

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

para_1_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_2_incorrect = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_control, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_1_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_1, 'w_s_all_recorded', w_actual_1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_2_correct = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_actual_2, 'w_s_all_recorded', w_actual_2, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial_1 = struct('R_t', 2);

para_Trivial_2 = struct('R_t', 5);

[~, ~, ~, ~, ~, Mean_0_default, Upper_0_default, Lower_0_default] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial_1);

[~, ~, ~, ~, ~, Mean_1_incorrect, Upper_1_incorrect, Lower_1_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_incorrect, para_Trivial_1);

[~, ~, I_1, ~, ~, Mean_1_correct, Upper_1_correct, Lower_1_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_1_correct, para_Trivial_1);

[~, ~, ~, ~, ~, Mean_2_incorrect, Upper_2_incorrect, Lower_2_incorrect] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_2_incorrect, para_Trivial_2);

[~, ~, I_2, ~, ~, Mean_2_correct, Upper_2_correct, Lower_2_correct] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para_2_correct, para_Trivial_2);

daysflip = [tau+1:total_time, total_time:-1:tau+1];


inBetween_0_default = [Lower_0_default(tau+1:total_time), fliplr(Upper_0_default(tau+1:total_time))];

inBetween_1_incorrect = [Lower_1_incorrect(tau+1:total_time), fliplr(Upper_1_incorrect(tau+1:total_time))];
inBetween_1_correct = [Lower_1_correct(tau+1:total_time), fliplr(Upper_1_correct(tau+1:total_time))];

inBetween_2_incorrect = [Lower_2_incorrect(tau+1:total_time), fliplr(Upper_2_incorrect(tau+1:total_time))];
inBetween_2_correct = [Lower_2_correct(tau+1:total_time), fliplr(Upper_2_correct(tau+1:total_time))];

inBetween_Area_1 = [Mean_1_correct(tau+1:total_time), fliplr(Mean_1_incorrect(tau+1:total_time))];

inBetween_Area_2 = [Mean_2_correct(tau+1:total_time), fliplr(Mean_2_incorrect(tau+1:total_time))];


%%
figure(1)

h(1) = plot(x, w_s_o, 'k');
hold on
h(2) = plot(x, w_s_a1, 'color', C(1, :));
h(3) = plot(x, w_s_a2, 'color', C(1, :), 'LineStyle', '--');

axis([0 15 0 0.45])

legend(h(1), 'Original')

title('Serial interval updates')
xlabel('Interval, $t$ (days)')
ylabel('Probability')

if Printer == 10

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_0.eps -eps

% 
% matlab2tikz('triagle.tex'); % black background

end
%%
figure(2)
clf
h(1) = plot([0 total_time], [para_Trivial_1.R_t para_Trivial_1.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
h(2) = fill(daysflip, inBetween_0_default, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_0_default(tau+1:total_time), 'color', 'k');

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 2]), '$R_t = 2$', 'Mean $\tilde{R}_t$', '95 \% CI')

title('Standard $R_t$ Inference')

if Printer == 10

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
export_fig Report_Constant_R_Disc_Serial_0_Default.eps -eps -opengl

end

%%
figure(3)
clf
h(1) = plot([0 total_time], [para_Trivial_1.R_t para_Trivial_1.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
% h(2) = fill(daysflip, inBetween_1_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);
% 
% h(3) = fill(daysflip, inBetween_1_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(2) = fill(daysflip, inBetween_Area_1', 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(4) = plot(tau+1:total_time, Mean_1_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_1_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_1_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2]), '$R_t = 2$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', 'Non-normalised cumulative error', 'Location', 'East')

title('Comparative inference when $R_t=2$ \& $\mu_{\mbox{\boldmath $\Delta$}}=0.7$')

if Printer == 2
figure(3)
% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [20 15]);
export_fig Report_Example_1_delta_eq_1.eps -eps -opengl

% 
% matlab2tikz('triagle.tex'); % black background

end

figure(4)
clf
h(1) = plot([0 total_time], [para_Trivial_2.R_t para_Trivial_2.R_t], 'color', [.5 .5 .5], 'LineWidth', 1.5);
hold on
% h(2) = fill(daysflip, inBetween_2_incorrect, C(4, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);
% 
% h(3) = fill(daysflip, inBetween_2_correct, C(1, :), 'LineStyle', 'none', 'FaceAlpha', 0.25);

h(2) = fill(daysflip, inBetween_Area_2, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25);


h(4) = plot(tau+1:total_time, Mean_2_incorrect(tau+1:total_time), 'color', C(4, :));
h(5) = plot(tau+1:switch_behaviour, Mean_2_correct(tau+1:switch_behaviour), 'color', C(1, :), 'LineStyle', '--');
h(6) = plot(switch_behaviour:total_time, Mean_2_correct(switch_behaviour:total_time), 'color', C(1, :));

h(7) = xline(switch_behaviour, '--', 'color', [.25 .25 .25]);

xlabel('Time, $t$ (days)')
ylabel('$\tilde{R}_t$ ')

legend(h([1 4 6 2]), '$R_t = 5$', 'Mean $\tilde{R}_t$ (original SI)', 'Mean $\tilde{R}_t$ (updated SI)', 'Non-normalised cumulative error', 'Location', 'East')

title('Comparative inference when $R_t=5$ \& $\mu_{\mbox{\boldmath $\Delta$}}=-0.7$')

if Printer == 2

% Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [20 15]);
export_fig Report_Example_2_delta_eq_1.eps -eps -opengl

end


