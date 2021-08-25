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

Raw_Data = readtable('../../../../Data/WHO-Ebola-CONFIDENTIAL/Line_List.csv');

patient_data = Raw_Data;

patient_data(:,1:5) = [];

Records = table2array(patient_data);

Times = table2array(Raw_Data(:, 3:5));

Active_Serial = zeros(1, size(Records, 1));

Delay = Serial;

for i = 1:size(Records,1)
    
    Delay(i) = datenum(Times(i, 3) - Times(i,1));
    
    if ismember(1, Records(i, :))
        
        Active_Serial(i) = find(Records(i,:) == 1, 1, 'first'); %Finds #days after tracking began that symptoms were suspected
        
    else
        
        Active_Serial(i) = find(Records(i,:) == 3, 1, 'first'); %Finds #days after tracking began that the case was confirmed
        
    end
    
end
   

Serial = Delay + Active_Serial;

figure(1)
clf
histogram(Serial, 'Normalization', 'probability')

ylabel('Probability')
xlabel('Delay between symptom onsets')

figure(2)
clf
histogram(Active_Serial, 'Normalization', 'probability')

ylabel('Probability')
xlabel('Delay between symptom onsets')
