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

N_o = 30;

pd = makedist('Normal','mu',18,'sigma',5);

trunc = truncate(pd, 0, inf);

x = linspace(1,N_o-1,N_o-1);

w_s_o = pdf(trunc,x)/sum(pdf(trunc,x));

w_s_o = [0 w_s_o];

% writematrix(w_s_o,'test_serial.csv') %For comparison to EpiEstim

w_s_a = w_s_o;

w_s_a(floor(0.4*N_o):end) = [];

w_s_a = w_s_a/sum(w_s_a);

N_a = length(w_s_a);

num_hybrids = max(N_a, N_o-N_a); %doesnt include the hybrid that is identical to original pmf

w_s_hy = cell(num_hybrids, 1); %Cell array since hybrids are different length

time_interest = [0 5 10 15 18];

for i = 1:num_hybrids  %doesnt includes the hybrid that is identical to old
    
    %If N_a is larger than N_o-i, then 
    
    length_hy = max(N_o-i, N_a);
    
    tmp = zeros(1, length_hy); %Length of hybrid
    
%     if N_o-i >= N_a+1 %The backtracking hasn't reached the new threshold
    
    tmp(1:end) = w_s_o(1:length_hy);
    
    if i <= N_a
        
        tmp(1:i) = w_s_a(1:i);
        
    elseif i > N_a
        
        tmp(1:N_a) = w_s_a(1:end);
        
    end
    
    
    
    tmp = tmp/(sum(tmp)); %Normalization
    
    w_s_hy{i} = tmp;
    
end

k = 0; %Used for pretty plots

for i = 1:num_hybrids
    figure(1)  
%     imshow(processo(:,:,1,i))
      j = i-1; %Used for easy indexing

      h(1) = plot(0:length(w_s_hy{i})-1, w_s_hy{i}, 'k');
      hold on
      h(2) = plot(0:length(w_s_a)-1, w_s_a, 'color', C(4, :), 'LineStyle', '--');
      h(3) = plot(0:length(w_s_o)-1, w_s_o, 'color', C(2, :), 'LineStyle', '--');  
        
      if N_o-i>=N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
        
          set(h(4),'linestyle','none')
          
          set(h(5),'linestyle','none')
          
          legend(h([3 2 1 4 5]), {'Original', 'Actual', 'Hybrid', 'Actual becomes relevant', 'Original still relevant'}, 'Location', 'best')
        
      elseif N_o-i>=N_a && i >N_a
        
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
          
          set(h(5),'linestyle','none')
          
          legend(h([3 2 1 5]), {'Original', 'Actual', 'Hybrid', 'Original still relevant'})
      
      elseif N_o-i<N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          set(h(4),'linestyle','none')
          
          legend(h([3 2 1 4]), {'Original', 'Actual', 'Hybrid', 'Actual becomes relevant'})
        
      end
      
      
      hold off
      
      
      title("Day "+i)
      
      axis([0 N_o 0 0.5])
      ylabel('Serial probability')
      xlabel('Interval, $t$ (Days)')
      
      
      F(i) = getframe(gcf) ;
      drawnow
      
    clf

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


%%



for i = 1:num_hybrids
    figure(1)  
%     imshow(processo(:,:,1,i))
      j = i-1; %Used for easy indexing

      h(1) = plot(0:length(w_s_hy{i})-1, w_s_hy{i}, 'k');
      hold on
      h(2) = plot(0:length(w_s_a)-1, w_s_a, 'color', C(4, :), 'LineStyle', '--');
      h(3) = plot(0:length(w_s_o)-1, w_s_o, 'color', C(2, :), 'LineStyle', '--');  
        
      if N_o-i>=N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
            
          set(h(4),'linestyle','none')
          
          set(h(5),'linestyle','none')
          
          legend(h([3 2 1 4 5]), {'Original', 'Actual', 'Hybrid', 'Actual becomes relevant', 'Original still relevant'}, 'Location', 'best')
        
      elseif N_o-i>=N_a && i >N_a
        
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
          
          set(h(5),'linestyle','none')
          
          legend(h([3 2 1 5]), {'Original', 'Actual', 'Hybrid', 'Original still relevant'})
      
      elseif N_o-i<N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          set(h(4),'linestyle','none')
          
          legend(h([3 2 1 4]), {'Original', 'Actual', 'Hybrid', 'Actual becomes relevant'})
        
      end
      
      
      hold off
      
      
      title("Day "+i)
      
      axis([0 N_o 0 0.5])
      ylabel('Serial probability')
      xlabel('Interval, $t$ (Days)')
      
      
      if ismember(j, time_interest)
          
          k = k+1;
          
          set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 15], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
          saveas(gcf, 'Hybrid_Example.eps')
          
          export_fig Hybrid_Example.eps -eps -r300 -painters -transparent
          
          num = num2str(k);
          
          movefile('Hybrid_Example.eps', [strcat('Hybrid_Example', num) '.eps'])
          
      end

end

%%

figure(2)
clf
for i = 1:num_hybrids
    
    j = i-1; %Used for easy indexing
    
    if ismember(j, time_interest)
      
      [~, q] = find(time_interest == j);
        
      subplot(1, 6, q)
      
      h(1) = plot(0:length(w_s_hy{i})-1, w_s_hy{i}, 'k');
      hold on
      h(2) = plot(0:length(w_s_a)-1, w_s_a, 'color', C(4, :), 'LineStyle', '--');
      h(3) = plot(0:length(w_s_o)-1, w_s_o, 'color', C(2, :), 'LineStyle', '--');  
        
      if N_o-i>=N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
            
          set(h(4),'linestyle','none')
          
          set(h(5),'linestyle','none')
        
      elseif N_o-i>=N_a && i >N_a
        
          h(5) = plot(N_o-i-1, w_s_hy{i}(end),'Marker', '.', 'MarkerSize', 25, 'color', C(2, :));
          
          set(h(5),'linestyle','none')
          
      elseif N_o-i<N_a && i <=N_a
        
          h(4) = plot(j, w_s_hy{i}(i),'Marker', '.', 'MarkerSize', 25, 'color', C(4, :));
          
          set(h(4),'linestyle','none')
          
      end
      
      title("Day "+i)
      
      if q == 3
      
          xlabel('Interval, $t$ (Days)')
      
      end
      
      if q == 1
         
          ylabel('Serial probability')
        
      end
        
    end
    


end

legend(h([3 2 1 4 5]), {'Original', 'Actual', 'Hybrid', 'Actual contribution', 'Original contribution'}, 'Position', [.82 .262 .1 .1])

set(gcf, 'Units', 'centimeters', 'Position', [0 0 30 8], 'PaperUnits', 'centimeters', 'PaperSize', [15 20]);
saveas(gcf, 'Hybrid_Progression.eps')
          
export_fig Hybrid_Progression.eps -eps -r300 -painters -transparent