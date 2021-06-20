
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

w_s_o = [0 0.1 0.2 0.3 0.2 0.1 0.05 0.03 0.02];

w_s_f = w_s_o;

% w_s_f(floor(0.7*N_o):end) = 0;
% 
% w_s_f = w_s_f/sum(w_s_f);

tau = 8;

para = struct('seed', 1996, 'total_time', 70, 'w_s_o', w_s_o, 'w_s_f', w_s_f, 'tau', tau, 'a', 1, 'b', 5, 'I_0', 10);

para_Trivial = struct('R_t', 2);

[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI_2('Fixed', 'Trivial', para, para_Trivial);

figure(1)
clf
plot(I)

figure(2)
clf
plot(Mean, 'r')
hold on
[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI_2('Perfect', 'Trivial', para, para_Trivial);
plot(Mean, 'b')

%% Now try for variable R_t

%Input a function

para_Variable = struct('R_t', @(x) 1.5 + sin(x));

[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI_2('Fixed', 'Variable', para, para_Variable);

figure(3)
clf
plot(I)

figure(4)
clf
plot(Mean, 'r')
hold on
[w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI_2('Perfect', 'Variable', para, para_Variable);
plot(Mean, 'b')
