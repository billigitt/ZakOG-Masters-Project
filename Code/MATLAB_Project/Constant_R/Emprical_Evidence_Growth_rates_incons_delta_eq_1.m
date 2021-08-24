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

Delta_Mat = reshape(Delta_Mat', [1 N 2*N+1])

Key = {'w_s_all_actual', 'R_t'};

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

R_start = 6;

R_end = 0.1;
 
days = 0:1:total_time(1);

para_Linear_Vary = struct('R_t', R_start + (R_end(1)-R_start(1))*days/total_time(1));

[Mean_Dif, Area_Dif, ~] = Sensitivity_Analysis_NonAbs(Key, para, para_Linear_Vary, 'Perfect', 'Trivial', 'Non-Hybrid');

[Abs_Mean_Dif, Abs_Area_Dif, ~] = Sensitivity_Analysis(Key, para, para_Linear_Vary, 'Perfect', 'Trivial', 'Non-Hybrid');
%% Plots
R_t = para_Linear_Vary.R_t;
figure(5)
clf

subplot(1, 2, 1)
imagesc(para_Linear_Vary.R_t, mu, Area_Dif)

xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{\int^{T_e}_{T_i} \tilde{R}_t^o(s) - \tilde{R}_t^a(s) \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}$', 'Interpreter', 'latex', 'FontSize', 18)

subplot(1, 2, 2)
imagesc(R_t, mu, Abs_Area_Dif)

xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}$', 'Interpreter', 'latex', 'FontSize', 18)

if Printer == 0

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 40 10], 'PaperUnits', 'centimeters', 'PaperSize', [40 10]);
saveas(gcf, 'Area_R_t_delta_eq_1.eps')
export_fig Area_R_t_delta_eq_1.eps -eps -r300 -painters -transparent

end



figure(6)
clf
subplot(1, 2, 1)
imagesc(R_t, mu, Mean_Dif)

xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{ \tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)}{\tilde{R}_t^a(T_e)}$', 'Interpreter', 'latex', 'FontSize', 18)


subplot(1, 2, 2)
imagesc(R_t, mu, Abs_Mean_Dif)


xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{ |\tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)|}{\tilde{R}_t^a(T_e)}$', 'Interpreter', 'latex', 'FontSize', 18)

if Printer == 0

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 40 10], 'PaperUnits', 'centimeters', 'PaperSize', [40 10]);
saveas(gcf, 'Mean_R_t_delta_eq_1.eps')
export_fig Mean_R_t_delta_eq_1.eps -eps -r300 -painters -transparent

end