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

%% Plots

t = linspace(-5, 15, 10000);

t_long = linspace(-7.5, 15, 10000);

options = optimset('TolFun',1e-40,'MaxFunEvals',1e8,'Maxiter',1e9);

epsilon = 1e-6;

fun1 = @(x) 1+sin(pi*(13*(exp(x+log(6/7)))/(6*(exp(x+log(6/7))+1))-1))  ;
Fun1 = @(x) sin(pi*(13*(exp(x+log(6/7)))/(6*(exp(x+log(6/7))+1))-1))  ;
fun2 = @(x) 1-sin(pi*(13*(exp(x+log(6/7)))/(6*(exp(x+log(6/7))+1))-1))  ;
Fun2 = @(x) -sin(pi*(13*(exp(x+log(6/7)))/(6*(exp(x+log(6/7))+1))-1))  ;
fun3 = @(x) 1+sin(pi*(13*(exp(x+log(12)))/(6*(exp(x+log(12))+1))-2))  ;
Fun3 = @(x) sin(pi*(13*(exp(x+log(12)))/(6*(exp(x+log(12))+1))-2))  ;
fun4 = @(x) 1-sin(pi*(13*(exp(x+log(12)))/(6*(exp(x+log(12))+1))-2))  ;
Fun4 = @(x) -sin(pi*(13*(exp(x+log(12)))/(6*(exp(x+log(12))+1))-2))  ;

fun5 = @(x) 1+sin(pi*(7*(exp(x+log(6)))/(6*(exp(x+log(6))+1))-1))  ;
fun6 = @(x) 1-sin(pi*(7*(exp(x+log(6)))/(6*(exp(x+log(6))+1))-1))  ;

fun7 = @(x) 1+cos(pi*(8*(exp(x+log(3)))/(3*(exp(x+log(3))+1))-2))  ;
fun8 = @(x) 1-cos(pi*(8*(exp(x+log(3)))/(3*(exp(x+log(3))+1))-2))

fun9 = @(x) -sin(pi*((exp(x))/((exp(x)+1))-0.5))  ;

r1 = fzero(Fun1, [epsilon 10], options);
r2 = fzero(Fun2, [epsilon 10], options);
r3 = fzero(Fun3, [-10 -epsilon], options);
r4 = fzero(Fun4, [-10 -epsilon], options);

options = optimset('TolFun',1e-10,'MaxFunEvals',1e8,'Maxiter',1e9);

y1 = zeros(1, length(t));
y2 = zeros(1, length(t));
y3 = zeros(1, length(t_long));
y4 = zeros(1, length(t_long));

y5 = zeros(1, length(t_long));
y6 = zeros(1, length(t_long));

y7 = zeros(1, length(t_long));
y8 = zeros(1, length(t_long));

y9 = zeros(1, length(t_long));

for i = 1:length(t)
    
   y1(i) = fun1(t(i));
   y2(i) = fun2(t(i));
    
end

for i = 1:length(t_long)
    
   y3(i) = fun3(t_long(i));
   y4(i) = fun4(t_long(i));
   y5(i) = fun5(t_long(i));
   y6(i) = fun6(t_long(i));
   
   y7(i) = fun7(t_long(i));
   y8(i) = fun8(t_long(i));
   
   y9(i) = fun9(t_long(i));
   
end

%% Plot 

f1 = figure(1)
clf
subplot(2, 2, 1)
plot(t, y1, 'k')

yline(0.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) t(end) -0.5 max(y1)+0.1])

plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)

plot(r1, 1, '^', 'MarkerFaceColor', C(4, :), 'MarkerSize', 10, 'color', C(4, :))

yticks([0.5 1])
set(gca,'YTickLabel', {'$w_1^a/w_1^o$','1'});
xticks([0])
ylabel('$h(r)$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -6;
hYLabel.Position(2) = 1.5;
title('$\mu_{\mbox{\boldmath $\Delta$}}<0$, $w_1^a< w_1^o \Rightarrow r^*>0$')
box on
subplot(2, 2, 3)

plot(t, y2, 'k')

yline(1.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) t(end) -0.5 max(y2)+0.1])

plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)

plot(r2, 1, '^', 'MarkerFaceColor', C(4, :), 'MarkerSize', 10, 'color', C(4, :))

yticks([1 1.5])
set(gca,'YTickLabel', {'1', '$w_1^a/w_1^o$'});
xticks([0])
ylabel('$h(r)$')
xlabel('$r$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -6;
hYLabel.Position(2) = 0.5;
hXLabel = get(gca,'XLabel');
hXLabel.Position(2) = -0.6;
title('$\mu_{\mbox{\boldmath $\Delta$}}>0$, $w_1^a> w_1^o \Rightarrow r^*>0$')
box on
subplot(2, 2, 4)

yticks([1 1.5])
set(gca,'YTickLabel', {'1', '$w_1^a/w_1^o$'});
xticks([0])

yline(1.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

plot(t_long, y3, 'k')

plot(r3, 1, '^', 'MarkerFaceColor', C(4, :), 'MarkerSize', 10, 'color', C(4, :))

plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)

axis([ t_long(1) t_long(end) -0.5 max(y3)+0.1])

xlabel('$r$')
hXLabel = get(gca,'XLabel');
hXLabel.Position(2) = -0.6;
title('$\mu_{\mbox{\boldmath $\Delta$}}<0$, $w_1^a> w_1^o \Rightarrow r^*<0$')
box on
subplot(2, 2, 2)

plot(t_long, y4, 'k')

yline(0.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([ t_long(1) t_long(end) -0.5 max(y4)+0.1])

h(1) = plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35);

h(2) = plot(r4, 1, '^', 'MarkerFaceColor', C(4, :), 'MarkerSize', 10, 'color', C(4, :))

xticks([0])
yticks([0.5 1])
set(gca,'YTickLabel', {'$w_1^a /w_1^o$','1'})

title('$\mu_{\mbox{\boldmath $\Delta$}}>0$, $w_1^a< w_1^o \Rightarrow r^*<0$')

legend(h([1 2]), '$h(0)=1$', '$h(r)=1$ $\mathrm{s.t.}$ $r \neq 0$')

% sgtitle('Qualitative behaviour when $\alpha = 2$ (two sign changes in \boldmath${\Delta}$)') 

ax11.FontSize = 18*20/40;
ax12.FontSize = 18*20/40;
ax13.FontSize = 18*20/40;
ax14.FontSize = 18*20/40;

box on

if Printer == 1

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 40 20], 'PaperUnits', 'centimeters', 'PaperSize', [40 20]);
saveas(gcf, 'Delta_Eq_2_Schematic.pdf')
% export_fig Delta_Eq_2_Schematic.eps -eps -r300 -painters -transparent

end

%% delta = 1

figure(2)
clf
ax21 = subplot(2, 1, 1);

yline(1.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

plot(t_long, y5, 'k')

axis([ t_long(1) t_long(end) -0.5 max(y5)+0.1])

xticks([0])
yticks([1 1.5])
set(gca,'YTickLabel', {'1', '$w_1^a /w_1^o$'})
ylabel('$h(r)$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
h(1) = plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35);
hYLabel.Position(2) = .5;
hYLabel.Position(1) = -9;
title('$\mu_{\mbox{\boldmath $\Delta$}}<0$')

box on
ax22 = subplot(2, 1, 2);

yline(0.5, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

hold on

yline(1, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

plot(t_long, y6, 'k')

axis([ t_long(1) t_long(end) -0.5 max(y6)+0.1])

xticks([0])
yticks([0.5 1])
set(gca,'YTickLabel', {'$w_1^a /w_1^o$','1'})

ylabel('$h(r)$')
xlabel('$r$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -9;
hYLabel.Position(2) = 1.5;
h(1) = plot(0, 1, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35);

hXLabel = get(gca,'XLabel');
set(hXLabel,'rotation',0,'VerticalAlignment','middle')
hXLabel.Position(2) = -0.75;

title('$\mu_{\mbox{\boldmath $\Delta$}}>0$')

% sgtitle('Qualitative behaviour when $\alpha = 1$ (one sign change in \boldmath${\Delta}$)') 

legend(h((1)), '$h(0)=1$', 'Location', 'NorthEast')

width=20;

ax21.FontSize = 18;
ax22.FontSize = 18;

box on

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 20], 'PaperUnits', 'centimeters', 'PaperSize', [20 20]);
saveas(gcf, 'Delta_Eq_1_Schematic.pdf')
% export_fig Delta_Eq_1_Schematic.eps -eps -r300 -painters -transparent

end
 %%
 
figure(3)
clf
ax1 = subplot(2, 1, 1)
plot(t_long, y7, 'k')

hold on

yline(2, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) 5 -0.5 max(y1)+.5])

h(1) = plot(0, 2, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)


yticks([0.5 2])
set(gca,'YTickLabel', {'$w_1^a/w_1^o$','1'});
xticks([0])
ylabel('$h(r)$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -5.5;
hYLabel.Position(2) = 1.25;
title("$\mu_{\mbox{\boldmath $\Delta$}}=0$, $\sigma^2_a< \sigma^2_o (\iff w_1^a<w_1^o)$:"+newline+"$\tilde{R}_t^o$ systematicly under-estimates $R_t$")
box on

ax2 = subplot(2, 1, 2)

plot(t_long, y8, 'k')
hold on
yline(0, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) 5 -0.5 max(y1)+0.5])

plot(0, 0, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)


yticks([0 1.5])
set(gca,'YTickLabel', {'1', '$w_1^a/w_1^o$'});
xticks([0])
ylabel('$h(r)$')
xlabel('$r$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -5.5;
hYLabel.Position(2) = .75;

hXLabel = get(gca,'XLabel');
% set(hXLabel,'rotation',0,'VerticalAlignment','middle')
hXLabel.Position(1) = -1.5;
hXLabel.Position(2) = -.65;

title("$\mu_{\mbox{\boldmath $\Delta$}}=0$, $\sigma^2_a> \sigma^2_o (\iff w_1^a>w_1^o)$:"+newline+"$\tilde{R}_t^o$ systematicly over-estimates $R_t$")

legend(h((1)), '$h(0)=1$', 'Location', 'best')

box on



ax1.FontSize = 18;
ax2.FontSize = 18;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 30], 'PaperUnits', 'centimeters', 'PaperSize', [20 30]);
saveas(gcf, 'Delta_Eq_2_Var_Sigma_Schematic.pdf')
% export_fig Delta_Eq_1_Schematic.eps -eps -r300 -painters -transparent

end

%%

figure(4)
clf
ax1 = axes('Position',[0.25 0.2 0.72 0.74]);
plot(t_long, -y9, 'k--')

hold on

plot(t_long, y9, 'k')

yline(0, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

% axis([t(1) 5 -0.5 max(y1)+.5])

% h(1) = plot(0, 2, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)


yticks([0])
set(gca,'YTickLabel', {'1'});
xticks([0])
set(gca,'XTickLabel', {'1'});

xlabel('True $R_t$')
set(gca, 'FontSize', 18)

ylabel('$\frac{\mathrm{Actual \: Inference}}{\mathrm{Ideal \: Inference}}$', 'FontSize', 22)


hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -4.12;
hYLabel.Position(2) = 0.4;

axis([-3 3 -1.1 1.1])

% title("$\mu_{\mbox{\boldmath $\Delta$}}=0$, $\sigma^2_a< \sigma^2_o (\iff w_1^a<w_1^o)$:"+newline+"$\tilde{R}_t^o$ systematicly under-estimates $R_t$")
box on





if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 18 12], 'PaperUnits', 'centimeters', 'PaperSize', [18 12]);
saveas(gcf, 'Wrong_Schematic_2.pdf')
% export_fig Delta_Eq_1_Schematic.eps -eps -r300 -painters -transparent

end


%%

figure(5)
clf
ax1 = subplot(1, 2, 1)
plot(t_long, y7, 'k')

hold on

yline(2, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) 5 -0.5 max(y1)+.5])

h(1) = plot(0, 2, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)

yticks([2])
set(gca,'YTickLabel', {'1'});
xticks([0])
set(gca,'XTickLabel', {'1'});
ylabel('$\frac{\tilde{R}^o_t}{\tilde{R}^a_t}$')
xlabel('True $R_t$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -5.5;
hYLabel.Position(2) = 1.25;
title("$\sigma^2_{\mathrm{true}}< \sigma^2_{\mathrm{recorded}}$"+newline+"$\tilde{R}_t^o$ systematicly under-estimates $R_t$")
box on

ax2 = subplot(1, 2, 2)

plot(t_long, y8, 'k')
hold on
yline(0, 'color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LineStyle', '--')

axis([t(1) 5 -0.5 max(y1)+0.5])

plot(0, 0, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 35)


yticks([0])
set(gca,'YTickLabel', {'1'});
xticks([0])
set(gca,'XTickLabel', {'1'});
ylabel('$\frac{\tilde{R}^o_t}{\tilde{R}^a_t}$')
xlabel('True $R_t$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -5.5;
hYLabel.Position(2) = 1.25;

hXLabel = get(gca,'XLabel');
set(hXLabel,'rotation',0,'VerticalAlignment','middle')
% hXLabel.Position(1) = -1.5;
hXLabel.Position(2) = -1.0;

title("$\sigma^2_{\mathrm{true}}> \sigma^2_{\mathrm{recorded}}$"+newline+"$\tilde{R}_t^o$ systematicly over-estimates $R_t$")

box on



ax1.FontSize = 18;
ax2.FontSize = 18;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 28 12], 'PaperUnits', 'centimeters', 'PaperSize', [28 12]);
saveas(gcf, 'Var_Sigma_Schem_Pres.pdf')
% export_fig Delta_Eq_1_Schematic.eps -eps -r300 -painters -transparent

end


%%
figure(6)
clf
axes('Position',[0.23 0.15 0.7 0.74]);

h = plot(t_long, y9+1, 'k--');

hold on

plot(t_long, y6, 'color', C(2, :))

plot(t_long, y8*0.5+1, 'color', C(2, :)*1/0.74)
plot(t, y2, 'color', C(2, :)*0.8)
axis([-5 5 -1.1 2.1])
yline(1, 'LineStyle', '--')
yticks([1])
set(gca,'YTickLabel', {'1'});
xticks([0])
set(gca,'XTickLabel', {'1'});


xlabel('True $R_t$')
set(gca, 'FontSize', 18)

ylabel('$\frac{\mathrm{Actual \: Inference}}{\mathrm{Ideal \: Inference}}$', 'FontSize', 22)

legend(h, "Current belief"+newline+"about estimation bias", 'Location', 'SouthWest')

hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -6.65;
hYLabel.Position(2) = 0.4;

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 18 12], 'PaperUnits', 'centimeters', 'PaperSize', [18 12]);
saveas(gcf, 'Summary_Schematic.pdf')
% export_fig Delta_Eq_1_Schematic.eps -eps -r300 -painters -transparent

end