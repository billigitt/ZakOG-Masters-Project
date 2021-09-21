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
%%

alph_1 = load('h_Linear_1.mat');

alph_2_Sigma_N_0 = load('h_Sigma_1.mat');

alph_2_Sigma_0 = load('h_Significant_1.mat');

mat_1 = alph_1.h_1;

mat_2_0 = alph_2_Sigma_0.Ratio;

mat_2_N_0 = alph_2_Sigma_N_0.Ratio;

R_start = 6;

R_end = 0.1;
 
days_p1 = 0:.1:100;

days_1  = 0:1:100;

R_t_1 = R_start + (R_end-R_start)*days_1/100;

R_t_p1 = R_start + (R_end-R_start)*days_p1/100;

figure(1)
clf
h(1) = plot(R_t_1, mat_1(1, :));

hold on

h(2) = plot(R_t_1, mat_1(end, :));

h(3) = plot(R_t_p1, mat_2_0(1, :));

h(4) = plot(R_t_p1, mat_2_0(end, :));

h(5) = plot(R_t_p1, mat_2_N_0(1, :));

h(6) = plot(R_t_p1, mat_2_N_0(end, :));

legend(h([1 2 3 4 5 6]), '1, 1', '1, end','a', 'b', 'c', 'd')