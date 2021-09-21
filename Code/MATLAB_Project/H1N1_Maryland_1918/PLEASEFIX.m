%Script for running and plotting SEIR ODE model 

%Cleaning
clc

close all

%Figure settings
set(0, 'defaultaxesfontsize', 16)
set(0, 'defaultlinelinewidth', 2)
set(0,'DefaultTextInterpreter', 'latex')
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
addpath('../Functions', '~/Documents/MATLAB/export_fig', '~/Documents/MATLAB/matlab2tikz')
savepath ../Functions/pathdef.m
%Defining model parameters, using day^-1 as units for beta, sigma, gamma and days for time
%NB: Recriprocal of 2 and 5 give correct parameters from data in question
para = struct('beta', 0.6,'sigma', 1/2, 'gamma', 1/5, 'N',8.982e6);

maxtime = 400;

%Define ICs
ICs = struct('S', para.N-1,'E', 0,'I', 1, 'R',0);

%Run model with no controls
[Classes] = ODE_SEIR_model_q1(para,ICs,maxtime);

%% Q1 (a): Plot dynamics

%Define colour spectrum
C  = [0.3686 0.3098 0.6353; 0.2005 0.5593 0.7380; 0.4558 0.7897 0.6458;...
    0.8525 0.2654 0.3082; 0.6196 0.0039 0.2588];

fn = fieldnames(Classes);


%% Q1 (b) (i): Compute final size of epidemic

%We round R_inf since  population is a whole number
R_inf = round(Classes.R(end));

%% Q1 (b) (iii): Compute expected outbreak duration

Epi_present = Classes.I>1;

% ix finds index of last infection that remains until Day 400 (Index 401)
ix = find(Epi_present,1, 'last');

% t_duration converts the index to time [index 1 = Day 0 etc.]
t_duration = ix-1;

%% Q1 (c): Plot of R_e(t)

R_e = para.beta*Classes.S/(para.N*para.gamma);

%Plot
f = figure(1);
clf



yyaxis left
h(1) = plot(Classes.t, R_e, 'color', 'k');
xlabel('Time, $t$ (days)','Interpreter','latex')
ylabel('True $R_t$','Interpreter','latex')
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% 
axis([0 120 0 3])
set(gca, 'FontSize', 20)
% 
% r = para.gamma*Classes.I.*(R_e-ones(length(R_e), 1));


yyaxis right
l(3) = plot(Classes.t, r, 'color', C(1, :))
xlabel('Time, $t$ (days)','Interpreter','latex')
ylabel("Epidemic growth"+newline+"rate, $r(t)$",'Interpreter','latex')

% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%  varies with time'},'Interpreter','latex')
% axis([0 120 -2e5 2.5e5])

ax = gca;
ax.YAxis(1).Color = 'k';

ax.YAxis(2).Color = C(1, :);


set(gca, 'FontSize', 20)





Printer=1;
if Printer == 1


set(gcf, 'Units', 'centimeters', 'Position', [0 0 22 11], 'PaperUnits', 'centimeters', 'PaperSize', [22 11]);
saveas(gcf, 'Intro_Conglom1.pdf')

end