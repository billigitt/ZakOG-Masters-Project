clc
clear all
close all

Key = {'tau', 'I_0'};

N_o = 15;

pd_o = makedist('Normal','mu',8,'sigma',2);
pd_a1 = makedist('Normal','mu',10,'sigma',2);

trunc_o = truncate(pd_o, 0, inf);
trunc_a1 = truncate(pd_a1, 0, inf);

x = linspace(1,N_o,N_o);

w_s_o = pdf(trunc_o,x)/sum(pdf(trunc_o,x));
w_s_a1 = pdf(trunc_a1,x)/sum(pdf(trunc_a1,x));

% Epidemiological and Inference parameters

total_time = 100;

tau = 2:6;

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

I_0 = 10:10:100;

para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_s_o, 'w_s_all_recorded', w_s_a1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);

R_start = 10;

R_end = 0.6;

days = 0:1:total_time(1);

para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));

[Mean_Dif, Area_Dif, ~] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Variable', 'Non-Hybrid')
