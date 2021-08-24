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
%% Start with analysis on how the serial inteval affects R_t inference

Key = {'tau', 'R_t'};

N_o = 10;

% pd_o = zeros(length(mu_vec));
x = linspace(1,N_o,N_o);

[~, d, ~] = Delta_Generator_Sigmoid_1(10, -.015, -0.3);

w_s_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

w_s_a1 = d + w_s_o;


% for i = 1:length(mu_vec)
% 
%     pd_o(i) = makedist('Normal','mu',mu_vec(i),'sigma',2);
%     
%     trunc_o(i) = truncate(pd_o(i), 0, inf);
% 
%     w_s_o(i, :) = pdf(trunc_o(i),x)/sum(pdf(trunc_o(i),x));
%         
% end

 %This is an important difference to old sensitivity analysis. We now only look at one

% Epidemiological and Inference parameters

total_time = 100;

tau = 3:10;

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

I_0 = 100;

para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', [w_s_o; w_s_a1], 'w_s_all_recorded', [w_s_o; w_s_a1], 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);

R_start = 4;

R_end = 0.9;

days = 0:1:total_time(1);

para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));

[Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Trivial', 'Non-Hybrid');

%%

figure(1)
clf
imagesc(para_Linear_Vary.R_t, tau, Area_Dif)
xlabel('True $R_t$')
ylabel('$\tau$')
% title('Mis-matching $\mathrm{Mean}(w^o)$ to $\mathrm{Mean}(w^a)$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s}{\int^{T_e}_{T_i} \tilde{R}_t^a(s) \mathrm{d}s}$', 'Interpreter', 'latex')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Area_Sensitivity_to_Rt_tau_Constant_R_t.eps')

export_fig Area_Sensitivity_to_Rt_tau_Constant_R_t.eps -eps -r300 -painters -transparent

end


figure(2)
clf
imagesc(para_Linear_Vary.R_t, tau, Mean_Dif)
xlabel('True $R_t$')
ylabel('$\tau$')
% title('Mis-matching $\mathrm{Mean}(w^o)$ to $\mathrm{Mean}(w^a)$')
set(gca,'YDir','normal')

c = colorbar;
colormap parula
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$\frac{|\tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)|}{\tilde{R}_t^a(T_e)}$', 'Interpreter', 'latex')

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Mean_Sensitivity_to_Rt_tau_Constant_R_t.eps')

export_fig Mean_Sensitivity_to_Rt_tau_Constant_R_t.eps -eps -r300 -painters -transparent

end


% %%
% 
% Key = {'w_s_all_actual', 'w_s_all_recorded'};
% 
% N_o = 15;
% 
% sigma_vec = 0.1:0.1:4;
% 
% % pd_o = zeros(length(mu_vec));
% x = linspace(1,N_o,N_o);
% 
% w_s_o = zeros(length(sigma_vec), N_o);
% 
% for i = 1:length(sigma_vec)
% 
%     pd_o(i) = makedist('Normal','mu',8,'sigma',sigma_vec(i));
%     
%     trunc_o(i) = truncate(pd_o(i), 0, inf);
% 
%     w_s_o(i, :) = pdf(trunc_o(i),x)/sum(pdf(trunc_o(i),x));
%         
% end
% 
% pd_a1 = pd_o;
% 
% trunc_a1 = trunc_o;
% 
% w_s_a1 = w_s_o;
% 
% % Epidemiological and Inference parameters
% 
% total_time = 100;
% 
% tau = 7;
% 
% switch_behaviour = 40;
% 
% delay = 0;
% 
% update_behaviour = switch_behaviour + delay;
% 
% I_0 = 100;
% 
% para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_s_o, 'w_s_all_recorded', w_s_a1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);
% 
% R_start = 5;
% 
% R_end = 0.6;
% 
% days = 0:1:total_time(1);
% 
% para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));
% 
% [Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Variable', 'Non-Hybrid');
% 
% %%
% 
% figure(3)
% clf
% imagesc(sigma_vec, sigma_vec, Area_Dif)
% xlabel('$\sigma_o$')
% ylabel('$\sigma_a$')
% title('Mis-matching $\mathrm{Var}(w^o)$ to $\mathrm{Var}(w^a)$')
% set(gca,'YDir','normal')
% % all_colors = get(0, 'DefaultAxesColorOrder');
% 
% % colormap(all_colors(1:2,:))
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Area_Sensitivity_to_Std_SI_Decreasing_R_t.eps')
% 
% export_fig Area_Sensitivity_to_Std_SI_Decreasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% 
% %%
% 
% figure(4)
% clf
% imagesc(sigma_vec, sigma_vec, Mean_Dif)
% xlabel('$\sigma_o$')
% ylabel('$\sigma_a$')
% title('Mis-matching $\mathrm{Var}(w^o)$ to $\mathrm{Var}(w^a)$')
% set(gca,'YDir','normal')
% 
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Mean_Sensitivity_to_Std_SI_Decreasing_R_t.eps')
% 
% export_fig Mean_Sensitivity_to_Std_SI_Decreasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% %%
% 
% Key = {'w_s_all_actual', 'w_s_all_recorded'};
% 
% N_o = 15;
% 
% mu_vec = 0:0.1:9;
% 
% % pd_o = zeros(length(mu_vec));
% x = linspace(1,N_o,N_o);
% 
% w_s_o = zeros(length(mu_vec), N_o);
% 
% for i = 1:length(mu_vec)
% 
%     pd_o(i) = makedist('Normal','mu',mu_vec(i),'sigma',2);
%     
%     trunc_o(i) = truncate(pd_o(i), 0, inf);
% 
%     w_s_o(i, :) = pdf(trunc_o(i),x)/sum(pdf(trunc_o(i),x));
%         
% end
% 
% pd_a1 = pd_o;
% 
% trunc_a1 = trunc_o;
% 
% w_s_a1 = w_s_o;
% 
% % Epidemiological and Inference parameters
% 
% total_time = 100;
% 
% tau = 7;
% 
% switch_behaviour = 40;
% 
% delay = 0;
% 
% update_behaviour = switch_behaviour + delay;
% 
% I_0 = 100;
% 
% para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_s_o, 'w_s_all_recorded', w_s_a1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);
% 
% R_start = 0.6;
% 
% R_end = 5;
% 
% days = 0:1:total_time(1);
% 
% para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));
% 
% [Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Variable', 'Non-Hybrid');
% 
% %%
% 
% figure(5)
% clf
% imagesc(mu_vec, mu_vec, Area_Dif)
% xlabel('$\mu_o$')
% ylabel('$\mu_a$')
% title('Mis-matching $\mathrm{Mean}(w^o)$ to $\mathrm{Mean}(w^a)$')
% set(gca,'YDir','normal')
% 
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Area_Sensitivity_to_Mean_SI_Increasing_R_t.eps')
% 
% export_fig Area_Sensitivity_to_Mean_SI_Increasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% 
% figure(6)
% clf
% imagesc(mu_vec, mu_vec, Mean_Dif)
% xlabel('$\mu_o$')
% ylabel('$\mu_a$')
% title('Mis-matching $\mathrm{Mean}(w^o)$ to $\mathrm{Mean}(w^a)$')
% set(gca,'YDir','normal')
% 
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Mean_Sensitivity_to_Mean_SI_Increasing_R_t.eps')
% 
% export_fig Mean_Sensitivity_to_Mean_SI_Increasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% 
% %%
% 
% Key = {'w_s_all_actual', 'w_s_all_recorded'};
% 
% N_o = 15;
% 
% sigma_vec = 0.1:0.1:4;
% 
% % pd_o = zeros(length(mu_vec));
% x = linspace(1,N_o,N_o);
% 
% w_s_o = zeros(length(sigma_vec), N_o);
% 
% for i = 1:length(sigma_vec)
% 
%     pd_o(i) = makedist('Normal','mu',8,'sigma',sigma_vec(i));
%     
%     trunc_o(i) = truncate(pd_o(i), 0, inf);
% 
%     w_s_o(i, :) = pdf(trunc_o(i),x)/sum(pdf(trunc_o(i),x));
%         
% end
% 
% pd_a1 = pd_o;
% 
% trunc_a1 = trunc_o;
% 
% w_s_a1 = w_s_o;
% 
% % Epidemiological and Inference parameters
% 
% total_time = 100;
% 
% tau = 7;
% 
% switch_behaviour = 40;
% 
% delay = 0;
% 
% update_behaviour = switch_behaviour + delay;
% 
% I_0 = 100;
% 
% para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_s_o, 'w_s_all_recorded', w_s_a1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);
% 
% R_start = 0.6;
% 
% R_end = 5;
% 
% days = 0:1:total_time(1);
% 
% para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));
% 
% [Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Variable', 'Non-Hybrid');
% 
% %%
% 
% figure(7)
% clf
% imagesc(sigma_vec, sigma_vec, Area_Dif)
% xlabel('$\sigma_o$')
% ylabel('$\sigma_a$')
% title('Mis-matching $\mathrm{Var}(w^o)$ to $\mathrm{Var}(w^a)$')
% set(gca,'YDir','normal')
% % all_colors = get(0, 'DefaultAxesColorOrder');
% 
% % colormap(all_colors(1:2,:))
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\int^{T_e}_{T_i} |\tilde{R}_t^o(s) - \tilde{R}_t^a(s)| \mathrm{d}s$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Area_Sensitivity_to_Std_SI_Increasing_R_t.eps')
% 
% export_fig Area_Sensitivity_to_Std_SI_Increasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% 
% %%
% 
% figure(8)
% clf
% imagesc(sigma_vec, sigma_vec, Mean_Dif)
% xlabel('$\sigma_o$')
% ylabel('$\sigma_a$')
% title('Mis-matching $\mathrm{Var}(w^o)$ to $\mathrm{Var}(w^a)$')
% set(gca,'YDir','normal')
% 
% c = colorbar;
% set(c,'TickLabelInterpreter','latex')
% ylabel(c, '$\tilde{R}_t^o(T_e) - \tilde{R}_t^a(T_e)$', 'Interpreter', 'latex')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'Mean_Sensitivity_to_Std_SI_Increasing_R_t.eps')
% 
% export_fig Mean_Sensitivity_to_Std_SI_Increasing_R_t.eps -eps -r300 -painters -transparent
% 
% end
% 
% 
% %%
% 
% figure(9)
% 
% perm = [1 3 2 4]
% ax = zeros(4,1);
% for i = 1:4
%     ax(i)=subplot(2,2,i);
% end
% % Now copy contents of each figure over to destination figure
% % Modify position of each axes as it is transferred
% for i = 1:4
%     figure(perm(i))
%     h = get(gcf,'Children')
%     if i ~=4
%         h(1) = [];
%     end
%     newh = copyobj(h,9)
%     for j = 1:length(newh)
% posnewh = get(newh(j),'Position');
% possub  = get(ax(i),'Position');
% set(newh(j),'Position',...
% [possub(1) possub(2) possub(3) possub(4)])
%     end
%     delete(ax(i));
% end
% 
% figure(9)