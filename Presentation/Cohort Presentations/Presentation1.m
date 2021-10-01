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
Printer = 0;
%%

% pd = makedist('Gamma','a',5,'b',1);
% 
% trunc = truncate(pd, 0, 7);
% 
x = 1:3;
y = 1:4;

w_s_o = [0.2 0.5 0.3];
% w_s_o_discrete = pdf(pd, y)/sum(pdf(pd, y));

% w_s_f = pdf(trunc, x);

I = [1 2 1 1];

I_Gen = 0.3*1+0.5*2+0.2*1;

%% Plots

hf = figure(1);
clf
% % scale_flip_w_s_o = 32/max(w_s_f)*fliplr(w_s_o);
% ax1 = subplot(2, 2, 1)
% bar([0 x], [0 w_s_o], 'BarWidth', 0.2, 'FaceColor', C(1, :))
% xticks([0 1 2 3])
% 
% 
% 
% % ha = annotation('arrow');
% % ha.Parent = hf.CurrentAxes;
% % ha.HeadStyle = 'deltoid';
% % ha.X = [x(45) 15];
% % ha.Y = [scale_flip_w_s_o(45) scale_flip_w_s_o(45)];
% % 
% % hb = annotation('arrow');
% % hb.Parent = hf.CurrentAxes;
% % hb.HeadStyle = 'deltoid';
% % hb.X = [15 x(45)];
% % hb.Y = [scale_flip_w_s_o(45) scale_flip_w_s_o(45)];
% 
% % set(gca, 'ytick', [])
% % xticks([0 x(45) 15])
% % xticklabels({'$T$', '$t$', '0'})
% 
% % ht = text(0.5*(15+x(45)), scale_flip_w_s_o(45)+5, {'Time between';'symptoms appearing';'in subsequent infections'}, 'interpreter', 'latex', 'FontSize', 13);
% % set(ht,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle')
% 
% ylabel("GI probability")
% 
% xlabel('Interval (days)')
% hYLabel = get(gca,'YLabel');
% hYLabel.Position(1) = -1.62;

% title("$R_t$ inference visualised via $\mbox{\boldmath $ I$}_t \cdot \mathrm{flip}(\mbox{\boldmath $ w$})$", 'Position', [4.6 0.65])

ax2 = subplot(2, 2, 2)

bar(y, [fliplr(w_s_o) 0], 'BarWidth', 0.2, 'FaceColor', C(1, :))

xticks([1 2 3 4])
xticklabels({'3', '2', '1', '0'})
xlabel('Interval (days)')

ylabel("GI probability")
hYLabel = get(gca,'YLabel');
hYLabel.Position(1) = -0.8;

ax = subplot(2, 2, 4)
ax.Clipping = 'off'
b1 = bar(y, I, 'FaceColor', [.75 .75 .75]);
hold on

% h=bar(y(4), I(4))
% set(h,'EdgeColor', 'red', 'FaceColor', [.5 .5 .5], 'LineWidth', 2);

b2 = bar(fliplr(w_s_o)*20/6, 'FaceColor', C(1, :), 'BarWidth', 0.2);

h=bar(1, fliplr(w_s_o(3))*20/6, 'FaceColor', C(1, :), 'BarWidth', 0.2)
% set(h,'EdgeColor', 'red', 'FaceColor', C(1, :), 'LineWidth', 2);

b3 = bar(4, I_Gen, 'BarWidth', 0.4, 'FaceColor', [0.9290 0.6940 0.1250]);
% h = bar(4, I_Gen, 'BarWidth', 0.4, 'FaceColor', [0.9290 0.6940 0.1250]);
% set(h,'EdgeColor', 'red', 'FaceColor', [0.9290 0.6940 0.1250], 'LineWidth', 2);


ylabel('Incidence')
xlabel('Time, $t$ (days)')

xticks([1 2 3 4])
xticklabels({'F', 'Sa', 'Su', 'M'})

legend([b1 b2 b3], {"Incidence", "Non-normalised"+newline+"GI", "$\mathrm{E}[I_t|R_t=1]$"}, 'Position', [0.2 0.62 0.15 0.25], 'FontSize', 18)



% ylabel("GI probability")
% 
% xlabel('Interval (days)')
% hYLabel = get(gca,'YLabel');
% hYLabel.Position(1) = -1.62;


ha = annotation('arrow');
ax.Clipping = 'off';


ha.Parent = hf.CurrentAxes;
ha.HeadStyle = 'plain';
ha.LineWidth = 10;
ha.HeadLength = 20;
ha.HeadWidth = 40;
% ha.X = [4.93 4.93];
% ha.Y = [2.8 2.2];
ha.Position = [5.23 2.8 0 -0.44]
ha.Color = [.5 .5 .5];

% hb = annotation('arrow');
% ax.Clipping = 'off';

% 
% hb.Parent = hf.CurrentAxes;
% hb.HeadStyle = 'plain';
% hb.LineWidth = 10;
% hb.HeadLength = 20;
% hb.HeadWidth = 40;
% % hb.X = [-1.6 0];
% % hb.Y = [4 4];
% hb.Position = [-1.8 4.05 0.85 0];
% hb.Color = [.5 .5 .5];
% set(hf, 'Position', [150 100 200 200])


ax1.FontSize = 18;
ax2.FontSize = 18;
ax.FontSize = 18;
Printer = 0;
if Printer == 0
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'R_t_Inference_Schematic_Norm.eps')

export_fig R_t_Inference_Schematic_Norm.eps -eps -r300 -painters -transparent

end

% figure(4)
% clf
% 
% I_end = zeros(length(I),1);
% 
% I_end(end) = I(end);
% 
% I_new = sum(I_Gen);
% 
% R_calc = zeros(length(I), 1);
% R_calc(end) = I_new;
% 
% ba = bar(y, [I' I_Gen'], 'stacked', 'FaceColor', 'flat');
% ba(1).CData = [0.75 0.75 0.75];
% ba(2).CData = C(4, :);
% hold on
% 
% ba(3) = bar(15, I(end), 'FaceColor', 'flat', 'LineWidth', 2, 'BarWidth', 0.9);
% 
% ba(3).CData = [0.5 0.5 0.5];
% 
% ba(4) = bar(15, I_new, 'FaceColor', 'flat', 'LineWidth', 2, 'BarWidth', 0.5);
% 
% ba(4).CData = C(4, :);
% 
% hold on
% 
% ba(5) = plot(x, 32/max(w_s_f)*fliplr([w_s_o]), 'color', C(1, :));
% 
% legend(ba([5 1 2 3 4]), {'Serial Interval*', 'Incidence', 'Contribution', 'Today''s Incidence', 'Projected Incidence'}, 'Location', 'Best')
% ylabel('Incidence')
% xlabel('Days')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'CohortPres1_4.eps')
% 
% export_fig CohortPres1_4.eps -eps -r300 -painters -transparent
% 
% end
% 
% figure(5)
% clf
% 
% h(1) = bar(y, I, 'FaceColor', [.75 .75 .75]);
% hold on
% 
% ix = find(w_s_f~=0, 1, 'last');
% 
% h(2) = plot(x, 32/max(w_s_f)*fliplr([w_s_o]), 'color', C(1, :));
% 
% h(3) = plot(x(end-ix+1:end), 32/max(w_s_f)*fliplr([w_s_f(1:ix)]), 'color', C(2, :));
% 
% legend(h([1 2 3]), {'Incidence', 'Original Serial*', 'Actual Serial*'}, 'Location', 'Best')
% ylabel('Incidence')
% xlabel('Days')
% %Save figure
% 
% if Printer == 1
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'CohortPres1_5.eps')
% 
% export_fig CohortPres1_5.eps -eps -r300 -painters -transparent
% 
% end
% 
% R_t_est = I(end)/I_new;
% 
