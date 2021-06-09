clc
clear all
close all



%TrivialEstimate_2

R_t = 3; %Fixed R_t

w_s = [0.1 0.2 0.3 0.2 0.1 0.05 0.03 0.02]; %Serial interval, e.g. odds
...of infecting after 1 day is 0.2.

total_time = 30; %Therefore total time+1 total pieces of data, since I_0
...is at time = 0.

days = 0:total_time;

tau = 5; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;

%Initial incidence

I = 1;

%Generate incidence data

I = 