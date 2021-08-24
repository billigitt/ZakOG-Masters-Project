clc
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

c = linspace(-1, 1, 100);

z = zeros(1, 100);

for i=1:100
z(i) = h_Transform(c(i), Delta_Matrix(1, :, 1));
end

plot(c, z)

fun = @(x) h_Transform(x, Delta_Matrix(1, :, 1));

options = optimset('TolFun',1e-40,'MaxFunEvals',1e8,'Maxiter',1e9, 'Display', 'iter');

c_star = fzero(fun,[0.00001 1], options)

%%
