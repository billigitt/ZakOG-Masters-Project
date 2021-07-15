%Cohort Presentation

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
Printer = 1;
%%

pd = makedist('Gamma','a',5,'b',1);

trunc = truncate(pd, 0, 7);

x = linspace(0,15,100);
y = linspace(0, 15, 16);

w_s_o = pdf(pd,x);
w_s_o_discrete = pdf(pd, y)/sum(pdf(pd, y));

w_s_f = pdf(trunc, x);

I = [1 1 3 8 7 10 12 13 11 16 22 19 21 24 28 32];

I_Gen = fliplr(w_s_o_discrete).*I;

%% Plots

hf = figure(1);
clf

scale_flip_w_s_o = 32/max(w_s_f)*fliplr(w_s_o);

plot(x, scale_flip_w_s_o, 'color', C(1, :))
hold on

ha = annotation('arrow');
ha.Parent = hf.CurrentAxes;
ha.HeadStyle = 'deltoid';
ha.X = [x(45) 15];
ha.Y = [scale_flip_w_s_o(45) scale_flip_w_s_o(45)];

hb = annotation('arrow');
hb.Parent = hf.CurrentAxes;
hb.HeadStyle = 'deltoid';
hb.X = [15 x(45)];
hb.Y = [scale_flip_w_s_o(45) scale_flip_w_s_o(45)];

set(gca, 'ytick', [])
xticks([0 x(45) 15])
xticklabels({'$T$', '$t$', '0'})

ht = text(0.5*(15+x(45)), scale_flip_w_s_o(45)+5, {'Time between';'symptoms appearing';'in subsequent infections'}, 'interpreter', 'latex', 'FontSize', 13);
set(ht,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle')

ylabel('Probability density')

xlabel('Days')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'CohortPres1_1.eps')

export_fig CohortPres1_1.eps -eps -r300 -painters -transparent

end

figure(2)
clf

bar(y, I, 'FaceColor', [.75 .75 .75])
xlabel('Days')
ylabel('Incidence')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'CohortPres1_2.eps')

export_fig CohortPres1_2.eps -eps -r300 -painters -transparent

end

figure(3)
clf

h(1) = bar(y, I, 'FaceColor', [.75 .75 .75]);
hold on

h(2) = plot(x, 32/max(w_s_f)*fliplr([w_s_o]), 'color', C(1, :));

ylabel('Incidence')
xlabel('Days')

legend(h([1 2]), {'Incidence', 'Serial Interval*'}, 'Location', 'Best')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'CohortPres1_3.eps')

export_fig CohortPres1_3.eps -eps -r300 -painters -transparent

end

figure(4)
clf

I_end = zeros(length(I),1);

I_end(end) = I(end);

I_new = sum(I_Gen);

R_calc = zeros(length(I), 1);
R_calc(end) = I_new;

ba = bar(y, [I' I_Gen'], 'stacked', 'FaceColor', 'flat');
ba(1).CData = [0.75 0.75 0.75];
ba(2).CData = C(4, :);
hold on

ba(3) = bar(15, I(end), 'FaceColor', 'flat', 'LineWidth', 2, 'BarWidth', 0.9);

ba(3).CData = [0.5 0.5 0.5];

ba(4) = bar(15, I_new, 'FaceColor', 'flat', 'LineWidth', 2, 'BarWidth', 0.5);

ba(4).CData = C(4, :);

hold on

ba(5) = plot(x, 32/max(w_s_f)*fliplr([w_s_o]), 'color', C(1, :));

legend(ba([5 1 2 3 4]), {'Serial Interval*', 'Incidence', 'Contribution', 'Today''s Incidence', 'Projected Incidence'}, 'Location', 'Best')
ylabel('Incidence')
xlabel('Days')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'CohortPres1_4.eps')

export_fig CohortPres1_4.eps -eps -r300 -painters -transparent

end

figure(5)
clf

h(1) = bar(y, I, 'FaceColor', [.75 .75 .75]);
hold on

ix = find(w_s_f~=0, 1, 'last');

h(2) = plot(x, 32/max(w_s_f)*fliplr([w_s_o]), 'color', C(1, :));

h(3) = plot(x(end-ix+1:end), 32/max(w_s_f)*fliplr([w_s_f(1:ix)]), 'color', C(2, :));

legend(h([1 2 3]), {'Incidence', 'Original Serial*', 'Actual Serial*'}, 'Location', 'Best')
ylabel('Incidence')
xlabel('Days')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'CohortPres1_5.eps')

export_fig CohortPres1_5.eps -eps -r300 -painters -transparent

end

R_t_est = I(end)/I_new;

