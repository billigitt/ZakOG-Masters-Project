
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

%% Data

incidence_data = readtable('../../../Data/WHO-COVID-19-global-data.csv');

idx_all = find(contains(incidence_data.Country, 'China')); %all indices for China

Dates = datestr(incidence_data.Date_reported(idx_all)); %All corresponding dates

Logical_start = Dates == '09-Jan-2020';
Logical_end = Dates == '13-Feb-2020';

idx_start = find(sum(Logical_start, 2) == length('09-Jan-2020'));
idx_end = find(sum(Logical_end, 2) == length('13-Feb-2020'));

idx = idx_all(idx_start:idx_end);

Incidence = incidence_data.New_cases(idx);

%% Serial Intervals

N_1 = 20; %Not sure if this is good

%Should think about whether these guys shrink or not.

N_2 = N_1; 

N_3 = N_1;

%

pd1 = makedist('Normal','mu',7.8,'sigma',5.2);

trunc1 = truncate(pd1, 0, inf);

x1 = linspace(1,N_1-1,N_1-1);

w_s_1 = pdf(trunc1,x1)/sum(pdf(trunc1,x1));

%

pd2 = makedist('Normal','mu',5.1,'sigma',5.0);

trunc2 = truncate(pd2, 0, inf);

x2 = linspace(1,N_2-1,N_2-1);

w_s_2 = pdf(trunc2,x2)/sum(pdf(trunc2,x2));

%

pd3 = makedist('Normal','mu',7.8,'sigma',4.6);

trunc3 = truncate(pd3, 0, inf);

x3 = linspace(1,N_3-1,N_3-1);

w_s_3 = pdf(trunc3,x3)/sum(pdf(trunc3,x3));

figure(1)

plot(x1, w_s_1)
xlabel('Interval, $t$ (Days)')
hold on
plot(x2, w_s_2)
plot(x3, w_s_3)

%% R_t Inference

tau = 7;

w_s_all_actual = [w_s_1; w_s_2; w_s_3];

w_s_all_recorded = [w_s_1; w_s_2; w_s_3];

switch_behaviour = [15 15+7];

delay = 0;

update_behaviour = switch_behaviour+ delay;

para = struct('seed', 196, 'total_time', length(Incidence), 'w_s_all_actual', w_s_all_actual, 'w_s_all_recorded', w_s_all_recorded, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', []);

para_Data = struct('Incidence', Incidence);
