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

%% Empirical evidence for how we can predict the non zero c value for which SI is not important

%Plan:

%Create a vector of mu_d (mean of Delta) running from -1 to 1. From this we
%find generate Delta-vectors. From that, we use fsolve to estimate c-values
%that yield h(c) = 1 where c =/= 0. Then we use the sensitivity analyses
%with key parameters, R and w_s_actual. We input w_s_actual with the
%(delta+w_o) values that we have already found (choosing a suitable w_o)
%and then (once we find a way to get R to correspond to c- should be
%doable?), we can then overlay the c-estimate with the sensitivity analyses
%and we should find they match up perfectly!

%% load('Large_Simulation_1') to skip running code- this will take roughly 8 hours!

Num_Serials = 100;

mu_d = linspace(.3, -.6, Num_Serials); %Differences in means of old and new serials

w_o = [0.0337 0.0387 0.0422 0.1437 0.1987 0.2040 0.1628 0.0739 0.0687 0.0336];

N = 10;

Days = 1:N;

Delta_Matrix = zeros(1, N, Num_Serials);

options = optimset('TolFun',1e-20,'MaxFunEvals',1e8,'Maxiter',1e9);

C_star = zeros(1, Num_Serials); %Vector of r-values that give h(r) = 1

R_star = zeros(1, Num_Serials); %Corresponding R-numbers (using eqn 3.6  in Wallinga-Lipsitch)

epsilon = 1e-10; %enables the domain to be roughly (0, x], i.e. non inclusive to 0

for i = 1:Num_Serials
    
    [Delta_Matrix(1, :, i), ~] = Delta_Generator_Quadratic_1(w_o, mu_d(i), +5); %Generates Deltas with different means. Set differing variance to 5 as enables good plot
    
    fun = @(x) h_Transform(x, Delta_Matrix(1, :, i))
    
    if fun(epsilon)<0 %There is a switch in the sign of mu_d which means the root will be found in a new place. Since fzero relies on a sign change we have to track that appropriatley
    
        C_star(i) = fzero(fun,[epsilon 10], options); %If we extend the range of mu_d, we may need larger than 10.
    
    else
    
        C_star(i) = fzero(fun,[-10 -epsilon], options);
        
    end
        
    R_star(i) = r_to_R(C_star(i), w_o + Delta_Matrix(1, :, i)); %Converts this c (r epidemic growth value) to an R number, using Wallinga-Lipsitch.
    
end

%% Sensitivity part

Key = {'w_s_all_actual', 'R_t'};

w_o = [0.0337 0.0387 0.0422 0.1437 0.1987 0.2040 0.1628 0.0739 0.0687 0.0336];

w_a_Matrix = repmat(w_o, 1, 1, Num_Serials)+ Delta_Matrix;

w_all_actual = zeros(2, N, Num_Serials);

w_all_actual(1, :, :) = repmat(w_o, 1, 1, Num_Serials);

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
 
days = 0:.1:total_time(1);

para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));

Num_Sims = 100;
    
[Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');

[Mean_Dif_NonAbs, Area_Dif_NonAbs, ~] = Sensitivity_Analysis_NonAbs(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');

%% Plots

R_t = para_Linear_Vary.R_t

%This figure shows the spread of Deltas that we will be studying. They are
%calculated based on what mu_d value we desire.

figure(1)
clf
plot(Days, Delta_Matrix(1, :, 1), 'color', C(1,:))

hold on

for i = 2:20
   
    h(1) = plot(Days, Delta_Matrix(1, :, i), 'color', C(1, :));
    
end

for i = 21:40
   
    h(2) = plot(Days, Delta_Matrix(1, :, i), 'color', C(2, :));
    
end

for i = 41:60
   
    h(3) = plot(Days, Delta_Matrix(1, :, i), 'color', C(3, :));
    
end

for i = 61:80
   
    h(4) = plot(Days, Delta_Matrix(1, :, i), 'color', C(4, :));
    
end

for i = 81:100
   
    h(5) = plot(Days, Delta_Matrix(1, :, i), 'color', C(5, :));
    
end

legend(h([1 2 3 4 5]),  strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu_d(1), '%4.2f'),'$,$', num2str(mu_d(20), '%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu_d(21), '%4.2f'),'$,$', num2str(mu_d(40), '%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu_d(41), '%4.2f'),'$,$', num2str(mu_d(60), '%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu_d(61), '%4.2f'),'$,$', num2str(mu_d(80), '%4.2f'), '$)$'), strcat('$\mu_{\mbox{\boldmath $\Delta$}} \in ($', num2str(mu_d(81), '%4.2f'),'$,$', num2str(mu_d(100), '%4.2f'), '$)$'), 'Location', 'best')

xlabel('Interval, $t$ (days)')
ylabel('$\mbox{\boldmath{$\Delta$}}(t)$')

hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -1.23;
title({'Distribution of';'quadratic \mbox{\boldmath$\Delta$} distributions for $\sigma^2_{\mbox{\boldmath$\Delta$}}$ \& $\alpha=2$'})

if Printer == 1

%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [20 15]);
saveas(gcf, 'Quadratic_Delta.eps')
export_fig Quadratic_Delta.eps -eps -r300 -painters -transparent

end

%This figure should tell us when to expect R>1 or R<1  and for what mu
%values. If the output is less than 0, it means that exactly one of 
%Delta_1 or mu_d are less than 0. Otherwise, they both have the same sign.
%The theory tells us that if they have the same sign, the output is
%positive and there will be R values >1 s.t h=1. If the output is negative,
%then the opposite is true and there are c/r values<0 (R<1) s.t h(c)=1.

figure(2)
clf
vec = zeros(1, Num_Serials);

for i = 1:Num_Serials

    vec(i)  = Delta_Matrix(1, 1, i);

end



plot(vec.*mu_d)

%This figure is less useful in the long run compared to R_star

figure(3)
clf
plot(mu_d, C_star)

figure(4)
clf
plot(mu_d, R_star)

f = figure(5);
clf
subplot(2, 11, [1 9])

imagesc(R_t, mu_d, Area_Dif.^(1/3))
colormap bone
colormap(f, flipud(colormap(f)))

% xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -.7;
set(gca,'YDir','normal')
title({'Inaccurate generation interval estimates may lead';'to accurate $R_t$ inference for non-trivial dynamics'})
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\sqrt{\frac{\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}}$', 'Interpreter', 'latex')

% legend(h, 'Numerical solve for $v(R(r)) \cdot  \mu_{\Delta}=0$', 'Location', 'SouthEast')

subplot(2, 11, [12 20])
imagesc(R_t, mu_d, Area_Dif.^(1/3))
colormap bone
colormap(f, flipud(colormap(f)))

hold on

h = plot(R_star, mu_d, 'k', 'LineStyle', '--');
xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -.7;

% hp4 = get(subplot(2,2,3),'Position');

c = colorbar('Position', [.8 .11 .05 0.816]);
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\bigg( \frac{\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s} \bigg)^{1/3}$', 'Interpreter', 'latex', 'FontSize', 18)

%legend(h, {"Numerical solve for"+newline+"$v(R(r)) \cdot  \mu_{\Delta}=0$"}, 'Location', 'SouthEast')


legend(h, {"Numerical solve for"+newline+" $v(R(r)) \cdot  \mu_{\mbox{\boldmath $\Delta$}}=0$"}, 'Location', 'SouthEast')


if Printer == 0
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 20], 'PaperUnits', 'centimeters', 'PaperSize', [20 20]);
saveas(gcf, 'Significant_Finding_1.eps')
export_fig Significant_Finding_1.eps -eps -r300 -painters -transparent

end

f = figure(6);
clf
% imagesc(R_t, mu_d, log2(1+tanh(10*Area_Dif_NonAbs)))

[X,Y] = meshgrid(R_t,mu_d);
H = surf(X, Y, log2(1+tanh(10*Area_Dif_NonAbs)))
colormap parula
colormap(f, flipud(colormap(f)))
H.EdgeColor = 'none';

hold on

plane = zeros(size(Area_Dif_NonAbs));

G = surf(X, Y, plane, 'FaceColor', 'k', 'FaceAlpha', 0.5);
G.EdgeColor = 'none';

% h = plot3(R_star, mu_d, 'k', 'LineStyle', '--');
xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\sqrt{\frac{ | \tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)|}{\tilde{R}_t^a(T_e)}}$', 'Interpreter', 'latex')

% legend(h, 'Numerical solve for $0 = \mu_{\Delta} \cdot $\boldmath$v$ $_R$', 'Location', 'SouthEast')

legend(h, {"Numerical solve for"+newline+"$v(R(r)) \cdot  \mu_{\Delta}=0$"}, 'Location', 'SouthEast')

% figure(10)
% 
% x = 0.0001-1:0.0001:10;
% 
% plot(x, log2(1+x), 'k')
% hold on
% plot(x, log(1+x), 'b')
% plot(x, log10(1+x), 'r')


