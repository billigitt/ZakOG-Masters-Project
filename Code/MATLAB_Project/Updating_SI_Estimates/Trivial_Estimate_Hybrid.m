
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

%% Run

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

w_s_a = w_s_o;

w_s_a(floor(0.7*N_o):end) = [];

w_s_a = w_s_a/sum(w_s_a);

N_a = length(w_s_a);

num_hybrids = max(N_a, N_o-N_a); %doesnt includes the hybrid that is identical to old

w_s_hy = cell(num_hybrids, 1);

for i = 1:num_hybrids  %doesnt includes the hybrid that is identical to old
    
    %If N_a is larger than N_o-i, then 
    
    length_hy = max(N_o-i, N_a);
    
    tmp = zeros(1, length_hy); %Length of hybrid
    
%     if N_o-i >= N_a+1 %The backtracking hasn't reached the new threshold
    
    tmp(1:end) = w_s_o(1:length_hy);
    
    tmp(1:i) = w_s_a(1:i);
    
    tmp = tmp/(sum(tmp)); %Normalization
    
    w_s_hy{i} = tmp;
    
end

for i = 1:num_hybrids
    figure(1)  
%     imshow(processo(:,:,1,i))
      
      plot(0:length(w_s_hy{i})-1, w_s_hy{i}, 'k')
      hold on
      plot(0:length(w_s_a)-1, w_s_a, 'color', C(4, :), 'LineStyle', '--')
      plot(0:length(w_s_o)-1, w_s_o, 'color', C(2, :), 'LineStyle', '--')
      hold off
      axis([0 10 0 0.5])
      F(i) = getframe(gcf) ;
      drawnow
end
  % create the video writer with 1 fps
  writerObj = VideoWriter('Hybrid_Progression.avi');
  writerObj.FrameRate = 1;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
