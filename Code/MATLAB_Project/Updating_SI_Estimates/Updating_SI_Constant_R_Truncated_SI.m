
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

w_s_f = w_s_o;

% w_s_f(floor(0.7*N_o):end) = 0;
% 
% w_s_f = w_s_f/sum(w_s_f);

tau = 7;

para = struct('seed', 1997, 'total_time', 100, 'w_s_o', w_s_o, 'w_s_f', w_s_f, 'tau', tau, 'a', 5, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI('Fixed', 'Trivial', para, para_Trivial);

figure(1)
clf
plot(I)
clear all
figure(2)
clf
plot(Mean, 'b')
hold on
[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI('Perfect', 'Trivial', para, para_Trivial);
plot(Mean, 'r')
