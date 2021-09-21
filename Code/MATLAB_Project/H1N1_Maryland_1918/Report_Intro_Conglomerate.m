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
% figure(1)
% for k = 1:4
% h(k) = plot(Classes.t, Classes.(fn{k}), 'color',C(k, :));
% hold on; 
% end
% 
% legend(h([1 2 3 4]), '$S(t)$', '$E(t)$', '$I(t)$', '$R(t)$','Interpreter','latex')
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('Frequency','Interpreter','latex')
% title({'A plot to show Infection'; 'dynamics with no interventions'},'Interpreter','latex')
% axis([0 400 0 1e7])

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

%Define R_e. This is identical to SIR model since the force of infection
%caused by the infected population is the same.
R_e = para.beta*Classes.S/(para.N*para.gamma);

%Plot
f = figure(1);
clf
subplot(3, 2, 1)
yyaxis left
h(1) = plot(Classes.t, R_e, 'color', 'k');
xlabel('Time, $t$ (days)','Interpreter','latex')
ylabel('$R_t$','Interpreter','latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel.Position(1) = -12;
% title({'A plot to show how the \emph{Instantaneous}'; '\emph{reproduction number} varies with time'},'Interpreter','latex')
axis([0 120 0 3])
hold on

%Define r
r = para.gamma*Classes.I.*(R_e-ones(length(R_e), 1));

%Plot
% figure(3)
yyaxis right
l(3) = plot(Classes.t, r, 'color', C(1, :))
xlabel('Time, $t$ (days)','Interpreter','latex')
ylabel('$r(t)$','Interpreter','latex')
hYLabel = get(gca,'YLabel');
hYLabel.Position(1) = 132;
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% title({'A plot of how the \emph{Infection}'; '\emph{growth gate} varies with time'},'Interpreter','latex')
axis([0 120 -2e5 2.5e5])

ax = gca;
ax.YAxis(1).Color = 'k';

ax.YAxis(2).Color = C(1, :);

% axes('position',[.15 .8 .1 .195])
% box on % put box around new pair of axes
% indexOfInterest = (Classes.t <= 75) & (Classes.t >= 74); % range of t near perturbation
% %plot(R_t(indexOfInterest),error_long_term_estimate(indexOfInterest), 'o', 'color', C(2, :)) % plot on new axes
% yyaxis left
% plot(Classes.t(indexOfInterest), R_e(indexOfInterest), 'color', 'k');
% hold on
% yyaxis right
% plot(Classes.t(indexOfInterest), r(indexOfInterest), 'color', C(1, :))
% ax = gca;
% ax.YAxis(1).Color = 'k';
% 
% ax.YAxis(2).Color = C(1, :);
% axis tight

%% Q2 
%Require change to ODE code , see ODE_SEIR_model_q2, to include new class and parameters

%Repeat parameterisation and ICs, introducing p_H, p_D and H class 
para2 = struct('beta', 0.6,'sigma', 1/2, 'gamma', 1/5, 'N',8.982e6, 'p_H', 0.01, 'p_D', 0.4);
ICs2 = struct('S', para.N-1,'E', 0,'I', 1, 'H', 0, 'R',0);
Classes2 = ODE_SEIR_model_q2(para2,ICs2,maxtime);

%% Q2 (b): Plot of hospitalisations and hospitalisation rate during epidemic

%Calculate daily_hospitalisation vector by staggering time and finding
%difference
daily_hospital = [Classes2.H;0]-[0;Classes2.H];
daily_hospital(end) = [];

%Plot
% figure(4)
% 
% %Subplot 1
% subplot(2, 1, 1 )
% plot(Classes2.t, Classes2.H, 'color', [0.5, 0.5, 0.5])
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel({'$H(t)$, Cumulative';'hospitalisations'},'Interpreter','latex')
% axis([0 t_duration 0 1e5])
% 
% %Subplot 2
% subplot(2, 1, 2 )
% plot(Classes2.t, daily_hospital, 'color', [0.5, 0.5, 0.5])
% ylabel('Daily hospitalisations','Interpreter','latex')
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% axis([0 t_duration 0 4000])

%% Q2 (c): Plot of deaths following hospitalisations

% Linear transformation of H to D
lagtime = zeros(14,1);
D = para2.p_D*[lagtime;Classes2.H];

% Extend time-doman to plot 2 weeks after 400 days
t_new = [Classes.t; (401:1:414)'];
daily_deaths = [D;0]-[0;D];
daily_deaths(end) = [];

%Plot
% figure(5)
% 
% %Subplot 1
% subplot(2, 1, 1)
% plot(t_new, D, 'color', C(5, :))
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel({'$D(t)$, Cumulative'; 'deaths'},'Interpreter','latex')
% axis([0 t_duration+14 0 4e4])
% 
% %Subplot 2
% subplot(2, 1, 2)
% plot(t_new, daily_deaths, 'color', C(5, :))
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('Daily deaths','Interpreter','latex')
% axis([0 t_duration+14 0 1600])

%% Q3
% Require change to ODE code, see ODE_SEIR_model_q3, so that mintime can be
% parsed
% Find first time I exceeds 100
% Add on a week and change parameter
% For pre-lockdown, ICs and parameters remain the same as Classes2 and ICs2

ICs_prelock = ICs2;
para_prelock = para2;

% Call daily hospitalisations vector (calculated in Q 2 (b)), find the time (index-1) of when
% hospitalisations first exceed 100. maxtime is 7 days later.

ix_infect100 = find(daily_hospital-100>0, 1, 'first');
time_infect100 = ix_infect100-1;

% Set parameters/ICs for pre lock-down (parameters same as normal)

mintime = 0;
maxtime = time_infect100+7;

Classes_prelock = ODE_SEIR_model_q3(para2,ICs2,mintime,maxtime);

% Set parameters/ICs for post lock-down

mintime = maxtime;
maxtime = mintime+51;
para_lock = struct('beta', 0.6*0.3,'sigma', 1/2, 'gamma', 1/5, 'N',8.982e6, 'p_H', 0.01, 'p_D', 0.4);
ICs_lock = struct('S', Classes_prelock.S(end),'E', Classes_prelock.E(end),'I', Classes_prelock.I(end), 'H', Classes_prelock.H(end), 'R',Classes_prelock.R(end));

Classes_lock = ODE_SEIR_model_q3(para_lock,ICs_lock,mintime,maxtime);

% Set parameters/ICs for post lock-down (parameters the same as pre lock-down)
mintime = maxtime;
%maxtime will now be the time that we calculated in 1 (b) (ii)
maxtime = t_duration;
para_postlock = para_prelock;
ICs_postlock = struct('S', Classes_lock.S(end),'E', Classes_lock.E(end),'I', Classes_lock.I(end), 'H', Classes_lock.H(end), 'R',Classes_lock.R(end));

Classes_postlock = ODE_SEIR_model_q3(para_postlock,ICs_postlock,mintime,maxtime);

%% Q 3 (a): Plot new infection dynamics from 0 to previous duration time

% %Plot
% figure(6)
% g(1) = plot(Classes_prelock.t, Classes_prelock.I, 'color', C(1, :));
% hold on
% g(2) = plot(Classes_lock.t, Classes_lock.I, 'color', C(2, :));
% hold on
% g(3) = plot(Classes_postlock.t, Classes_postlock.I, 'color', C(3, :));
% title({'Infection dynamics during'; 'different social restrictions'}, 'Interpreter', 'latex')
% legend(g([1 2 3]), 'Pre-Lockdown', 'Lockdown', 'Post-Lockdown', 'Interpreter', 'latex', 'Location', 'NorthWest')
% axis([0 181 0 1.1e6])
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('$I(t)$, Infected population','Interpreter','latex')

%% Q 3 (b): Plot of how R_e is impacted by the restrictions

%Define R_e for restrictions
R_e_prelock = para_prelock.beta*Classes_prelock.S/(para_prelock.N*para_prelock.gamma);
R_e_lock = para_lock.beta*Classes_lock.S/(para_lock.N*para_lock.gamma);
R_e_postlock = para_postlock.beta*Classes_postlock.S/(para_postlock.N*para_postlock.gamma);

t1 = Classes_prelock.t;
t1(end) = [];
t2 = Classes_lock.t;
t2(end) = [];

%Plot
% figure(7)
subplot(3, 2, 2)

l(1) = plot(Classes.t, R_e, 'k');
hold on
l(2) = plot(Classes_prelock.t, R_e_prelock, 'color', C(1, :));
hold on
plot(Classes_lock.t, R_e_lock, 'color', C(1, :))
hold on
plot(Classes_postlock.t, R_e_postlock, 'color', C(1, :))
xline(t1(end)+1, '--')
hold on
xline(t2(end)+1, '--')
hold on
% text1 = text(t1(end)-4,1, 'PHM begins', 'Interpreter', 'latex', 'FontSize', 18);
% hold on
% text2 = text(t2(end)-4,1, 'PHM ends', 'Interpreter', 'latex', 'FontSize', 18);
% set(text1,'Rotation',90)
% set(text2,'Rotation',90)
xlabel('Time, $t$ (days)','Interpreter','latex')
ylabel('$R_t$','Interpreter','latex')
hYLabel = get(gca,'YLabel');
hYLabel.Position(1) = -10;
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% title({'A plot to show how $R_e(t)$ is';'affected by social restrictions'},'Interpreter','latex')
legend(l([1 2]), '$R_t$ without PHM','$R_t$ with PHM', 'Interpreter','latex', 'Location', 'best')
axis([0 100 0 3])

% %% Q 3 (c): Final size of outbreak, R_inf & duration of outbreak
% 
% %We first need to extend our simulation time so that the infection has been
% %removed from the population, then we re-run the simulation
% maxtime = 600;
% Classes_postlocklong = ODE_SEIR_model_q3(para_postlock,ICs_postlock,mintime,maxtime);
% 
% %To calculate the final size of the outbreak, we need to know the sum of
% %the people recovered and hospitalised at the final time
% outbreaksize = round(Classes_postlocklong.R(end)+Classes_postlocklong.H(end));
% 
% %Find index of first 0 that remains until Day 600
% Epi_present_restrictions = Classes_postlocklong.I>1;
% ix_postlock = find(Epi_present_restrictions,1, 'last')+1;
% 
% %ix_restrictions counts the index of the whole simulation up to when I=0
% %The -2 is included so that we don't double count
% ix_restrictions = length(Classes_prelock.t)+length(Classes_lock.t)-2+ix_postlock;
% 
% % We then convert to days
% t_duration_restrictions = ix_restrictions-1;
% 
% %% Q 3(d): Daily incidence of hospitalisations and deaths
% 
% %Remove overlap between different periods
% t1 = Classes_prelock.t;
% t1(end) = [];
% H1 = Classes_prelock.H;
% H1(end) = [];
% t2 = Classes_lock.t;
% t2(end) = [];
% H2 = Classes_lock.H;
% H2(end) = [];
% t3 = Classes_postlocklong.t;
% t3(end) = [];
% H3 = Classes_postlocklong.H;
% H3(end) = [];
% 
% %Concatenate vectors to make plotting easier
% t_restrictions = [t1;t2;t3];
% H_restrictions = [H1; H2; H3];
% D_restrictions = para_prelock.p_D*[lagtime; H_restrictions];
% D_restrictions(end-13:end) = [];
% 
% %Calculate daily hospitalisations & deaths with restrictions (analagous to Q2 (b))
% daily_hospital_restrictions = [H_restrictions;0]-[0;H_restrictions];
% daily_hospital_restrictions(end) = [];
% 
% daily_deaths_restrictions = [D_restrictions;0]-[0;D_restrictions];
% daily_deaths_restrictions(end) = [];
% 
% %Calculate # Deaths
% D_donothing = para2.p_D*[lagtime; Classes2.H];
% total_deaths_donothing = round(D_donothing(end));
% total_deaths_restrictions = round(D_restrictions(end));
% 
% %Plot
% new  = figure(8);
% f(1) = plot(t_restrictions, daily_deaths_restrictions, 'color', C(5, :));
% hold on
% f(2) = plot(t_restrictions, daily_hospital_restrictions, 'color', [0.5, 0.5, 0.5]);
% hold on
% xline(t1(end)+1, '--')
% hold on
% xline(t2(end)+1, '--')
% hold on
% f(3) = plot(250, 5000, 'w');
% 
% axis([0 268 0 4000])
% text1 = text(t1(end)-4,1000, 'Lockdown begins', 'Interpreter', 'latex', 'FontSize', 18);
% text2 = text(t2(end)-4,1000, 'Lockdown ends', 'Interpreter', 'latex', 'FontSize', 18);
% set(text1,'Rotation',90)
% set(text2,'Rotation',90)
% 
% %Delete cells to make vector lengths match
% t_new(end-13:end)=[];
% daily_deaths(end-13:end) = [];
% 
% hold on
% f(4) = plot(t_new, daily_deaths, 'color', C(5, :), 'LineStyle', ':');
% hold on
% f(5) = plot(t_new, daily_hospital, 'color', [0.5, 0.5, 0.5], 'LineStyle', ':');
% hold on
% f(6) = plot(250, 5000, 'w');
% 
% lgd1 = legend(f([6 5 4 3 2 1]), '\textbf{Without restrictions}', 'Daily hospitalisations','Daily deaths in hospitals','\textbf{With restrictions}', 'Daily hospitalisations','Daily deaths in hospitals', 'Interpreter', 'latex');
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('Frequency','Interpreter','latex')
% title({'A plot to show how daily deaths \&'; 'hospitalisations vary during social restrictions'},'Interpreter', 'latex')
% 
% %% Q 4(a)
% 
% %Find maximum value of I, and corresponding time
% [max_I, ix_max] = max(Classes.I);
% t_peak_I = ix_max-1;
% 
% %Calculate herd immunity by substituting time found into R(t)
% herd_immunity = Classes.R(t_peak_I+1);
% 
% %Calculate overshoot
% overshoot = round(R_inf-herd_immunity);
% 
% %Plot
% figure(9)
% for k = 3:4
% h(k) = plot(Classes.t, Classes.(fn{k}), 'color',C(k, :));
% hold on; 
% end
% xline(t_peak_I, '--')
% hold on
% yline(herd_immunity, '--')
% text1 = text(t_peak_I-3,5.5e6, 'Peak infection time', 'Interpreter', 'latex', 'FontSize', 12);
% text2 = text(10, herd_immunity+2e5, 'Herd Immunity threshold', 'Interpreter', 'latex', 'FontSize', 12);
% text3 = text(120,6e6, 'Overshoot', 'Interpreter', 'latex', 'FontSize', 14);
% set(text1,'Rotation',90)
% 
% legend(h([3 4]), '$I(t)$', '$R(t)$','Interpreter','latex')
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('Frequency','Interpreter','latex')
% title({'A plot demonstrating the difference'; 'between herd immunity and outbreak size'},'Interpreter','latex')
% axis([0 t_duration 0 1e7])
% 
% %Arrow
% an = annotation('doublearrow');
% an.Parent = gca;
% an.Position  = [150, herd_immunity+1e5, 0, R_inf-herd_immunity-2e5];
% 
% %% Q 4(b)
% 
% %Set up parameters, ICs and time interval with new X parameter, percentage
% %of population with initial immunity
% para4 = struct('beta', 0.6,'sigma', 1/2, 'gamma', 1/5, 'N',8.982e6, 'p_H', 0.01, 'p_D', 0.4, 'X', 0.4);
% ICs4 = struct('S',para4.N*(1-para4.X)-1, 'E', 0, 'I', 1, 'H', 0, 'R', para4.N*para4.X);
% maxtime = 400;
% 
% Classes4 = ODE_SEIR_model_q2(para4,ICs4,maxtime);
% 
% %Outbreak size is sum of Recovered and Hospitalised communities
% outbreak_with_immunity = round(Classes4.R(end)+Classes4.H(end));
% 
% %Plot
% figure(10)
% p(1) = plot(Classes4.t, Classes4.I, 'color', [0.5, 0.5, 0.5]);
% hold on
% p(2) = plot(Classes.t, Classes2.I, 'k');
% 
% xlabel('Time, $t$ (\emph{days})','Interpreter','latex')
% ylabel('$I(t)$, Infected Population','Interpreter','latex')
% legend(p([2 1]), 'No initial immunity', '40\% initial immunity','Interpreter','latex')
% title({'Comparison of infection dynamics with'; 'and without some initial immunity'}, 'Interpreter', 'latex')
% axis([0 300 0 2e6])

%% Run H1N1_Maryland_1918

%Read data
set(0,'DefaultTextInterpreter', 'latex')
inference_data = readtable('Incidence_H1N1_Maryland_1918.csv');
serial_data = readtable('Serial_Interval_H1N1_Maryland_1918.csv');

local = inference_data.local;
imported = inference_data.imported;
total_time = length(inference_data.local)-1; % -1 because we define the start from

I = imported'+local'; %Prelimary calculation. NEED TO CHANGE.
%time 0

days = 0:total_time;

w_s = (serial_data.x)'; %Serial interval, e.g. odds
...of infecting after 1 days is 1/3.

tau = 7; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;



%Generate incidence data

Shape = zeros(1, total_time+1);
Scale = zeros(1, total_time+1);
Mean = zeros(1, total_time+1);
Upper = zeros(1, total_time+1);
Lower = zeros(1, total_time+1);

for t = tau+1:total_time+1 %tau+1 (Day tau) is the first you can start 
    
    ...from. t is the index, NOT the day.
%     disp(I(t-tau:t))
    Shape(t) = a + sum(I(t-tau+1:t)); % IS +1 RIGHT?? shape param at time t = i
    Scale(t) = 0;
    
    %Calculate summation of Lambdas
    
    for k  = t-tau+1:t %Support material suggests that we only look at 
        ...these values, contrary to Thompson-Stockwin paper
        
        I_relevant = I(1:k); %definitely want from k down, may want all
        ...the way down to 1
        
        Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, w_s); 
        %Takes I_k, I_{k-1} down to length of Serial
        %The 0 is included again. Previously it was ommitted, since then we
        %were *generating* the *following* day so we wanted the 'current'
        %(but really previous day) to be included. Here we are actually
        %thinking about the current day, i.e. k and so we re-include the 0,
        %to illustrate that the serial must be 0 at this value.
    end
    
%     I(1) = I_0; %First case is imported
    
    Scale(t) = 1/(Scale(t)+(1/b));
    
    Mean(t) = Scale(t)*Shape(t); %Mean of gamma distn is product of two 
    ...pars
    
    Upper(t) = gaminv(0.975, Shape(t), Scale(t));
    
    Lower(t) = gaminv(0.025, Shape(t), Scale(t));
    
end

%% Plots



%Plot Mean
ax = subplot(3, 2, [3 4])
% daysflip = [tau:total_time, total_time:-1:tau];
% inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
% h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1,'edgecolor', [1 1 1]);
% 
% hold on
% 
% h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k');
% 
% h(3) = plot([tau total_time], [1 1], 'k--', 'LineWidth', 1);
% 
% %Labels
% 
% % plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
% % h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% % % plot([days(1) days(end)], [a*b a*b], 'b--')
% 
% title({'Real-time $R_t$ inference of 1918 H1N1';['outbreak in Maryland, USA ($\tau =$ ', num2str(tau), ' days$)$']})
% ylabel('$\bar{R}_t$')
% xlabel('Time, $t$ (days)')
% 
% %Legend
% 
% legend(h([1 2]), '$\bar{R}_t$', '95\% Posterior CI', 'Location', 'Best')
% 
% Printer = 0;

% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'H1N1_Maryland_1918_Estimate_CI.eps')
% 
% export_fig H1N1_Maryland_1918_Estimate_CI.eps -eps -r300 -painters -transparent
% 
% end

% figure(2)
% clf

daysflip = [tau:total_time, total_time:-1:tau];
inBetween = [Lower(tau+1:total_time+1), fliplr(Upper(tau+1:total_time+1))];
yyaxis left

h(4) = bar(days, I, 'FaceColor', C(5, :), 'LineStyle', 'none');

ylabel('Incidence')

yyaxis right

h(1) = fill(daysflip, inBetween, [.75 .75 .75], 'LineWidth', 0.1,'edgecolor', [1 1 1]);

hold on

h(2) = plot(tau:total_time, Mean(tau+1:total_time+1), 'k', 'LineStyle', '-');

h(3) = plot([tau total_time], [1 1], 'k--', 'LineWidth', 1);

% title({'Real-time $R_t$ inference of 1918 H1N1';['outbreak in Maryland, USA ($\tau =$ ', num2str(tau), ' days$)$']})
ylabel('$\tilde{R}_t$')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

xlabel('Time, $t$ (days)')

ax = gca;
ax.YAxis(1).Color = C(5, :);

ax.YAxis(2).Color = 'k';

%Labels

% plot([days(1) days(end)], [gaminv(0.975, a, b) gaminv(0.975, a, b)], 'color', C(2, :), 'LineStyle', '--')
% h(4) = plot([days(1) days(end)], [gaminv(0.025, a, b) gaminv(0.025, a, b)], 'color', C(2, :), 'LineStyle', '--');
% % plot([days(1) days(end)], [a*b a*b], 'b--')



%Legend

legend(h([2 1 4]), '$\tilde{R}_t$', '95\% Posterior CI', 'Incidence Data', 'Position', [0.5999 0.525 0.1 0.105], 'Interpreter', 'latex')
% axis([0 total_time+tau 0 max(Upper)])

% if Printer == 0
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'H1N1_Maryland_1918_Estimate_CI_Zoom.eps')
% 
% export_fig H1N1_Maryland_1918_Estimate_CI_Zoom.eps -eps -r300 -opengl
% 
% end
ax.Position = [ 0.13 0.4 0.757 0.25];
%%

% figure(3)
% clf
% h(1) = bar(days, I, 'FaceColor', [.5 .5 .5]);
% % hold on
% % plot(days, R_t.^days, 'r')
% title('Incidence data for H1N1 outbreak in Maryland, USA')
% ylabel('Incidence')
% xlabel('Time, $t$ (days)')
% 
% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'H1N1_Maryland_1918_Estimate_Incidence.eps')
% 
% export_fig H1N1_Maryland_1918_Estimate_Incidence.eps -eps -r300 -painters -transparent
% 
% end
% 
% ix = find(w_s~=0, 1, 'last');
% 
% W_s = w_s(1:ix);
% 
% figure(4)
% plot(0:length(W_s)-1, W_s, 'k')
% title('EpiEstim Serial Interval')
% ylabel('Probability')
% xlabel('Time, $t$ (days)')

% if Printer == 1
% %Save figure
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
% saveas(gcf, 'H1N1_Maryland_1918_Estimate_Certain_Serial.eps')
% 
% export_fig H1N1_Maryland_1918_Estimate_Certain_Serial.eps -eps -r300 -painters -transparent
% 
% end

%% Comparison to EpiEstim

% Go to https://shiny.dide.imperial.ac.uk/epiestim/ and upload data as well
% as downloading EstimatedR.csv to compare to shape and scale parameters
% through time. The curves should lie exactly on top of one another.

inference_data = readtable('Estimated_R_H1N1_Maryland_1918.csv');

mu = inference_data.Mean_R_;

sigma = inference_data.Std_R_;

k = (mu.^2./sigma.^2);
theta = mu./k;
% x = 0:0.001:6;
% figure(5)
% clf
% plot(x, gampdf(x, k(end), theta(end)))

Printer = 1;


% figure(5)
% clf
subplot(3, 2, 5)
clf
h(1) = plot(k,'color', C(2, :));
hold on
h(2) = plot(Shape(tau+1:end), 'k', 'LineStyle', '--');
% title('Shape parameter comparison')
ylabel('$k$')
xlabel('Time, $t$ (days)')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
yticks([0 500 1500 2000])
subplot(3, 2, 6)
h(1) = plot(theta, 'color', C(2, :));
hold on
h(2) = plot(Scale(tau+1:end), 'k', 'LineStyle', '--');
% title('Scale parameter comparison')
ylabel('$\theta$')
xlabel('Time, $t$ (days)')
legend(h([1 2]), 'EpiEstim inference', 'Our inference', 'Position', [0.6305 0.24 0.05 0.07], 'Interpreter', 'latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
yticks([0 0.01 0.03 0.04])
% title('Typical $R_t$ and $r(t)$ behaviour, $\tilde{R}_t$ inference and Epi-Estim comparisons', 'Position', [.5 .5 .5])

if Printer == 1
%Save figure
set(gcf, 'Units', 'centimeters', 'Position', [60 10 40 20], 'PaperUnits', 'centimeters', 'PaperSize', [40 20]);
saveas(gcf, 'Intro_Conglom.eps')

% export_fig Intro_Conglom.eps -eps -r300 -painters -transparent

end
% f.Position = [2000 100 1040 800];