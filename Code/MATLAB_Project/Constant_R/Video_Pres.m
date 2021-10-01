%% Cleaning
clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 25)
set(0, 'defaultlinelinewidth', 3)
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
C  = [0.3686 0.3098 0.6353; 0.2005 0.5593 0.7380; 0.4558 0.7897 0.6458;...
    0.8525 0.2654 0.3082; 0.6196 0.0039 0.2588];
addpath('../Functions', '~/Documents/MATLAB/export_fig', '~/Documents/MATLAB/matlab2tikz', '~/Documents/MATLAB/ffmpeg')
savepath ../Functions/pathdef.m
Printer = 1;
%% Video_Pres


alph_2_Sigma_0 = load('h_Significant_1.mat');

mat_1 = alph_2_Sigma_0.Ratio;

Deltas = load('Deltas.mat');

Delta_Matrix = Deltas.Delta_Matrix;

R_start = 6;

R_end = 0.1;

days_1  = 0:.1:100;

R_t_1 = R_start + (R_end-R_start)*days_1/100;

figure('units','pixels','position',[0 0 1920 1080]) 
% For me it was a resolution issue. I build the figures at higher resolution and it works
subplot(1, 2, 1)
xlabel('Interval, $t$ (days)')
ylabel('Mis-specification, $\mbox{\boldmath $\Delta$}$')
ax1 = gca();
subplot(1, 2, 2)
xlabel('True $R_t$')
ylabel('Over/under-estimation $\tilde{R}_t$ bias')
ax2 = gca();
v = VideoWriter('Msc_Pres.avi');
v.FrameRate = 5;
v.Quality = 100;
open(v)

for t=100:-1:1
    if t == 100
      subplot(1, 2, 1)
      
      h = plot(ax1, 1:10, Delta_Matrix(1, :, 1), 'color', C(5, :));
      set(ax1, 'XLimMode', 'manual', 'YLimMode', 'manual');
      xlabel('Interval, $t$ (days)')
      ylabel('Mis-specification, $\mbox{\boldmath $\Delta$}$')
      axis([0 10 -0.13 0.21])
      subplot(1, 2, 2)
      g = plot(ax2, R_t_1, mat_1(1, :), 'color', C(2, :));
      set(ax2, 'XLimMode', 'manual', 'YLimMode', 'manual');
      xlabel('True $R_t$')
      ylabel('$\tilde{R}_t$ estimation bias')
      axis([0 2 0.983 1.05])
      yline(1, 'LineWidth', 2, 'LineStyle', '--')
    else
        subplot(1, 2, 1)
      set(h, 'YData', Delta_Matrix(1, :, t), 'color', C(5, :));
      axis([0 10 -0.13 0.21])
      subplot(1, 2, 2)
      set(g, 'YData', mat_1(t, :), 'color', C(2, :));
      axis([0 2 0.983 1.05])
      yline(1, 'LineWidth', 2, 'LineStyle', '--')
    end
    drawnow();
    writeVideo(v,getframe(gcf))
  end
  close(gcf)
  close(v)
  
  pathVideoMP4 = regexprep('Msc_Pres.avi','\.avi','.mp4'); % generate mp4 filename
  
   [~,~] = system(sprintf('ffmpeg -i %s -y -an -c:v libx264 -crf 0 -preset slow %s','Msc_Pres.avi','Msc_Pres.mp4')); % for this to work, you should have installed ffmpeg and have it available on PATH