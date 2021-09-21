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

%% Misspecification

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

figure(1)
clf
ax1 = axes('Position',[0.25 0.2 0.72 0.74]);
h(1) = plot(w_s_o, 'color', C(1, :), 'LineStyle', '-', 'Marker', '*', 'MarkerSize', 10);
hold on

h(2) = plot(w_s_a1, 'color', C(1, :), 'LineStyle', '--', 'Marker', '.', 'MarkerSize', 25);

h(3) = plot(w_s_a1-w_s_o, 'color', C(4, :), 'Marker', '.', 'MarkerSize', 25);

yline(0, 'LineWidth', 1.5)

legend(h([1 2 3]), 'Estimated GI', 'True GI', 'Mis-specification, $  \mbox{\boldmath $\Delta$}$', 'Location', 'SouthEast')

set(gca, 'FontSize', 18)

yticks([-0.1 0 0.1 0.2])

ylabel('Probability')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -2

xlabel('Interval, $t$ (days)')
Printer=1
if Printer == 1

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 28 15], 'PaperUnits', 'centimeters', 'PaperSize', [28 15]);
saveas(gcf, 'Misspecification.eps')
export_fig Misspecification.eps -eps -r300 -painters -transparent

end