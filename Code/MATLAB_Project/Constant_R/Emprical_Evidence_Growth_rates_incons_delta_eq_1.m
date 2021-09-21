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

M=210;

B = [N+1 2; 2*N+1 3];

A = blkdiag(B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B,B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B,B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B,B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B,B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B,B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B, B);

mu = linspace(1, -1, M);

x = zeros(2*M, 1);

for i = 1:M
   
    x(2*i-1) = 0;
    
    x(2*i) = 6*mu(i)/(N*(N+1));
    
end

phi = A\x;

a = phi(1:2:2*M-1);

b = phi(2:2:2*M);

Delta_Mat = repmat(b, 1, N) + a.*(1:N);

Delta_Mat = reshape(Delta_Mat', [1 N M]);

Key = {'w_s_all_actual', 'R_t'};

w_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

w_a_Matrix = repmat(w_o, 1, 1, M)+ Delta_Mat;

w_all_actual = zeros(2, N, M);

w_all_actual(1, :, :) = repmat(w_o, 1, 1, M);

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

[Mean_Dif, Area_Dif, ~] = Sensitivity_Analysis_NonAbs(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');

[Abs_Mean_Dif, Abs_Area_Dif, ~] = Sensitivity_Analysis(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');

[h_1, ~] = Sensitivity_Analysis_H(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');
%% Plots
R_t = para_Linear_Vary.R_t;
f = figure(1);
% clf
% 
% subplot(1, 2, 1)
% imagesc(para_Linear_Vary.R_t, mu, Area_Dif)
% 
% xlabel('True $R_t$')
% ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
% set(gca,'YDir','normal')
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% colormap parula
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\frac{\int^{T_e}_{T_i} \tilde{R}_t^o(s) - \tilde{R}_t^a(s) \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}$', 'Interpreter', 'latex', 'FontSize', 18)
% 
% subplot(1, 2, 2)
imagesc(R_t, mu, (Abs_Area_Dif).^(1/3))

xlabel('True $R_t$')
ylabel('Mean misspecification, $\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')
set(gca,'YDir','normal')
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
colormap bone
colormap(f, flipud(colormap(f)))
d = colorbar;
set(d,'TickLabelInterpreter','latex')
ylabel(d, '$\bigg(\frac{\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}\bigg)^{1/3}$', 'Interpreter', 'latex', 'FontSize', 14)
Printer=1;
set(gca, 'FontSize', 18)
Printer=1
% title('$R_t=1 \iff$ zero $\:$ $\:$ bias in $R_t$ estimation, for $\alpha=1$')
if Printer == 1

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20.5 12], 'PaperUnits', 'centimeters', 'PaperSize', [20.5 12]);
saveas(gcf, 'Area_R_t_delta_eq_1.eps')
export_fig Area_R_t_delta_eq_1.eps -eps -r300 -painters -transparent



end
Printer=1;

C  = [0.3686 0.3098 0.6353; 0.2005 0.5593 0.7380; 0.4558 0.7897 0.6458;...
    0.8525 0.2654 0.3082; 0.6196 0.0039 0.2588];
%%
f = figure(4);
clf
% imagesc(R_t, mu_d, log2(1+tanh(10*Area_Dif_NonAbs)))

[X,Y] = meshgrid(R_t,mu);
[c, h] = contourf(X, Y, h_1,[0.7 0.8 0.9 0.97 1 1.03 1.1 1.2 1.3], 'ShowText', 'on', 'LineColor', 'none');
set(h,'LineColor','flat')
clabel(c,h,'Color','k', 'Interpreter', 'latex');
colormap parula
% colormap(f, flipud(colormap(f)))
% H.EdgeColor = 'none';

% hold on

% plane = ones(size(h_1));
% 
% G = surf(X, Y, plane, 'FaceColor', 'k', 'FaceAlpha', 0.5);
% G.EdgeColor = 'none';

% h = plot3(R_star, mu_d, 'k', 'LineStyle', '--');
xlabel('True $R_t$')
ylabel('Mean misspecification, $\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

% hXLabel = get(gca,'YLabel');
% % hXLabel.Position(1) = 0;
% % hXLabel.Position(2) = 1;
% set(hXLabel,'rotation',0,'VerticalAlignment','middle')

% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$h(r; \mbox{\boldmath $\Delta$})$', 'Interpreter', 'latex')

%title('Shape of $h(R_t)$ matches with mathematical analysis, for $\alpha=1$')

legend(h, '$h(R_t)$')

set(gca, 'FontSize', 18)
Printer=0
if Printer == 0

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [20 15]);
saveas(gcf, 'h_delta_eq_1.eps')
export_fig h_delta_eq_1.eps -eps -r300 -painters -transparent

end
Printer=0;
%%
figure(2)
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
Printer=1
if Printer == 1

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 40 10], 'PaperUnits', 'centimeters', 'PaperSize', [40 10]);
saveas(gcf, 'Mean_R_t_delta_eq_1.eps')
export_fig Mean_R_t_delta_eq_1.eps -eps -r300 -painters -transparent

end

Days = 1:N;

figure(3)
clf
h(1) = plot(Days, Delta_Mat(1, :, 1), 'color', .5*C(1, :));

hold on

for i = 2:70
   
    plot(Days, Delta_Mat(1, :, i), 'color', .5*C(1, :))
    
end

for i = 71:140
   
    h(2) = plot(Days, Delta_Mat(1, :, i), 'color', 1.1*C(1, :));
    
end

for i = 141:210
   
    h(3) = plot(Days, Delta_Mat(1, :, i), 'color', C(1, :)*1/0.65);
    
end

legend(h([1 2 3]),  strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu(1),'%4.2f'),'$,$', num2str(mu(70),'%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu(71),'%4.2f'),'$,$', num2str(mu(140),'%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu(141),'%4.2f'),'$,$', num2str(mu(210),'%4.2f'), '$)$'), 'Location', 'North')

xlabel('Interval, $t$ (days)')
ylabel('Mis-specification, $\mbox{\boldmath{$\Delta$}}_t$')

% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% title({'Distribution of';'linear \mbox{\boldmath$\Delta$} distributions for $\alpha=1$'})

set(gca, 'FontSize', 18)
Printer=0
if Printer == 0

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [20 15]);
saveas(gcf, 'Linear_Delta.eps')
export_fig Linear_Delta.eps -eps -r300 -painters -transparent

end