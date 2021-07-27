
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
%% Run Updating_SI_Constant_R_Truncated_SI

%We create a make-believe disease with a serial interval that suddenly
%truncates. We assume the initial serial interval is a Normal+ distribution
%with mean 5 and std 5.

N_o = 12;

pd = makedist('Normal','mu',5,'sigma',2);

trunc = truncate(pd, 0, inf);

x = linspace(1,N_o-1,N_o-1);

w_s_o = pdf(trunc,x)/sum(pdf(trunc,x));

w_s_o = [0 w_s_o];

% writematrix(w_s_o,'test_serial.csv') %For comparison to EpiEstim

w_s_f = w_s_o;

w_s_f(floor(0.7*N_o):end) = 0;

w_s_f = w_s_f/sum(w_s_f);

tau = 8;

para = struct('seed', 1996, 'total_time', 70, 'w_s_o', w_s_o, 'w_s_f', w_s_f, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_cont_update_SI('Fixed', 'Trivial', para, para_Trivial);

% writematrix(I','test_data.csv') %For comparison to EpiEstim

figure(1)
clf
plot(I)

figure(2)
clf
h(1) = plot(Mean, 'r');
hold on
[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_cont_update_SI('Perfect', 'Trivial', para, para_Trivial);
h(2) = plot(Mean, 'b');

legend(h([1 2]), 'Fixed', 'Perfect', 'Location', 'Best')

%% Comparison to EpiEstim [CONFIRMED]

% Go to https://shiny.dide.imperial.ac.uk/epiestim/ and upload data as well
% as downloading EstimatedR.csv to compare to shape and scale parameters
% through time. The curves should lie exactly on top of one another.

% data = readtable('tester.csv');
% 
% mu = data.Mean_R_;
% 
% sigma = data.Std_R_;
% 
% k = (mu.^2./sigma.^2);
% theta = mu./k;
% x = 0:0.001:6;
% 
% figure(3)
% clf
% h(1) = plot(k, 'g');
% hold on
% h(2) = plot(Shape(tau+1:end), 'b--');
% title('Shape parameter comparison')
% ylabel('$k$')
% xlabel('Time, $t$ (days)')
% legend(h([1 2]), 'EpiEstim', 'My inference', 'Location', 'NorthWest')
% 
% 
% 
% figure(4)
% clf
% h(1) = plot(theta, 'g');
% hold on
% h(2) = plot(Scale(tau+1:end), 'r--');
% title('Scale parameter comparison')
% ylabel('$\theta$')
% xlabel('Time, $t$ (days)')
% legend(h([1 2]), 'EpiEstim', 'My inference')



%% Now try for variable R_t

%Input a function
para = struct('seed', 196, 'total_time', 365*7, 'w_s_o', w_s_o, 'w_s_f', w_s_f, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);
para_Variable = struct('R_t', @(x) 2.0 + 1.5*sin(2*pi*x/365));

[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_cont_update_SI('Fixed', 'Variable', para, para_Variable);

figure(3)
clf
plot(I)

figure(4)
clf
plot(Mean, 'r')
hold on
[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_cont_update_SI('Perfect', 'Variable', para, para_Variable);
plot(Mean, 'b')

y = 1:para.total_time;

plot(y, 2.0 + sin(2*pi*y/365))

%% Discrete changes- 1 quarantining

%The corresponding figure demonstrates that we can have quite strange
%behaviour for R_t inference when the serial intervals are not alligned
%with actual and recorded. We may even be surprised to see that these poor
%estimates are still present for a short period of time even with exact
%alignment between actual and recorded SIs (set delay=0). This happens
%because to switch between SIs instantly foes not take into account the
%fact that there will still be some infection flowing in from *previous*
%behaviours *after* an intervention.

w_s_all_actual = [w_s_o; w_s_f; w_s_o; w_s_f; w_s_o; w_s_f; w_s_o];

w_s_all_recorded = [w_s_o; w_s_f; w_s_o; w_s_f; w_s_o; w_s_f; w_s_o];

switch_behaviour = [30 50 70 110 160 165];

delay = 0;

update_behaviour = switch_behaviour+ delay;

para = struct('seed', 196, 'total_time', 300, 'w_s_all_actual', w_s_all_actual, 'w_s_all_recorded', w_s_all_recorded, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial = struct('R_t', 2);

[w_s_actual, w_s_recorded, I, Shape, Scale, Mean, Upper, Lower] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para, para_Trivial);

%We want to compare results with the hybrid SI method. To start with, we
%will just use one switch time.

switch_behaviour_hy = 50;
update_behaviuour_hy = 50;


%%
%%Plot


x1 = [0 0 30 30];
y1 = [0 2.5 2.5 0];

x2 = [60 60 70 70];
y2 = [0 2.5 2.5 0];

x3 = [120 120 160 160];
y3 = [0 2.5 2.5 0];

x4 = [165 165 170 170];
y4 = [0 2.5 2.5 0];

x5 = [175 175 300 300];
y5 = [0 2.5 2.5 0];


x7 = [30 30 40 40];
y7 = [0 2.5 2.5 0];

x8 = [70 70 80 80];
y8 = [0 2.5 2.5 0];

x9 = [160 160 165 165];
y9 = [0 2.5 2.5 0];



x10 = [40 40 50 50];
y10 = [0 2.5 2.5 0];

x11 = [80 80 110 110];
y11 = [0 2.5 2.5 0];



x6 = [170 170 175 175];
y6 = [0 2.5 2.5 0];

x12 = [50 50 60 60];
y12 = [0 2.5 2.5 0];

x13 = [110 110 120 120];
y13 = [0 2.5 2.5 0];
figure(5)
clf
% fill(x1,y1,'r')
% hold on
% fill(x2,y2,'r')
% fill(x3,y3,'r')
% fill(x4,y4,'r')
% fill(x5,y5,'r')
% fill(x6,y6,'y')
% fill(x7,y7,'b')
% fill(x8,y8,'b')
% fill(x9,y9,'b')
% fill(x10,y10,'g')
% fill(x11,y11,'g')
% fill(x12,y12,'y')
% fill(x13,y13,'y')

plot(Mean, 'k')



title('Mean when there are updates on SIs')



figure(6)
clf
plot(I)
title('Incidence')

%% Plot with Hybrids!

%We must try a much simpler approach (to start off with)

w_s_all_actual = [w_s_o; w_s_f];

w_s_all_recorded = [w_s_o; w_s_f];

switch_behaviour = 50;

delay = 0;

update_behaviour = switch_behaviour+ delay;

para = struct('seed', 196, 'total_time', 300, 'w_s_all_actual', w_s_all_actual, 'w_s_all_recorded', w_s_all_recorded, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial = struct('R_t', 2);

[w_s_actual_Hybrid, w_s_recorded_Hybrid, ~, ~, ~, Mean_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Hybrid', para, para_Trivial);

[w_s_actual_Non_Hybrid, w_s_recorded_Non_Hybrid, ~, ~, ~, Mean_Non_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para, para_Trivial);


figure(7)

h(1) = plot(Mean_Non_Hybrid, 'k');
hold on
h(2) = plot(Mean_Hybrid, 'r');

legend(h([1 2]), {'No Hybrid SIs', 'Hybrid SIs'}, 'Location', 'best')

ylabel('$R_t$ inference')
xlabel('Time, $t$ (Days)')

%%

w_s_all_actual(2, :) = w_s_all_actual(1, :);

w_s_all_recorded(2, :) = w_s_all_recorded(1, :);

switch_behaviour = 50;

delay = 0;

update_behaviour = switch_behaviour+ delay;

para = struct('seed', 196, 'total_time', 300, 'w_s_all_actual', w_s_all_actual, 'w_s_all_recorded', w_s_all_recorded, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 1000);

para_Trivial = struct('R_t', 2);

[w_s_actual_Hybrid, w_s_recorded_Hybrid, ~, ~, ~, Mean_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Hybrid', para, para_Trivial);

[w_s_actual_Non_Hybrid, w_s_recorded_Non_Hyrbid, ~, ~, ~, Mean_Non_Hybrid, ~, ~] = R_infer_disc_update_SI('Perfect', 'Trivial','Non-Hybrid', para, para_Trivial);

clf
figure(8)
plot(Mean_Hybrid, 'k')
hold on
plot(Mean_Non_Hybrid, 'r--')
