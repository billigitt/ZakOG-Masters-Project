clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 2)
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Set2'))

%%

%TrivialEstimate

R_t = 2; %Fixed R_t

w_s = [0.1 0.2 0.3 0.2 0.1 0.05 0.03 0.02]; %Serial interval, e.g. odds
...of infecting after 1 day is 0.2.

% w_s = [1/8 0 2/8 1/8 1/8 1/8 1/8 1/8];
% 
% w_s = [1/8 0 2/8 5/8];
% 
% w_s = [0.25 0.25 0.25 0.25];
% 
% w_s = [0 0 1 0 0]

total_time = 50; %Therefore total time+1 total pieces of data, since I_0
...is at time = 0.

days = 0:total_time;

tau = 5; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;

%Initial incidence

I = 1;

%Generate incidence data

for t = 1:total_time
   
    I_new = poissrnd(R_t*Incidence_Generator(I, w_s));
    I = [I_new, I];
    
end

I = fliplr(I); %Now it is old to new, i.e. I_0...I_total_time

%Calulcate shape and scale posterior parameters for Gamma distn

shape = zeros(1, total_time+1);
scale = zeros(1, total_time+1);
mean = zeros(1, total_time+1);
upper = zeros(1, total_time+1);
lower = zeros(1, total_time+1);

for t = tau+1:total_time+1 %tau+1 (Day tau) is the first you can start 
    ...from. t is the index, NOT the day.
   
    shape(t) = a + sum(I(t-tau:t)); %shape param at time t = i
    scale(t) = 0;
    
    for k  = t-tau:t %Index t-tau to t, i.e. day t-tau-1 to t-1
        
        I_switch = I(1:k); %definitely want from k down
        
        I_switch = fliplr(I_switch); %I_{k-1},...I_0
        
        scale(t) = scale(t) + Incidence_Generator(I_switch, w_s); %Takes I_k, 
        ...I_{k-1} down to length of Serial

    end
    
    scale(t) = 1/(scale(t)+(1/b));
    
    mean(t) = scale(t)*shape(t); %Mean of gamma distn is product of two 
    ...pars
    
    upper(t) = gaminv(0.995, shape(t), scale(t));
    
    lower(t) = gaminv(0.005, shape(t), scale(t));
    
end

figure(1)
% plot(days, upper)
% hold on
% plot(days, lower)
plot(days, mean)

hold on
daysflip = [days, fliplr(days)];
inBetween = [lower, fliplr(upper)];
I_switch = fill(daysflip, inBetween, 'r', 'FaceAlpha', 0.25, 'edgealpha', 0);
title({['Trivial $R_t$ (=', num2str(R_t), ') estimate with'];['$w_s = [$',num2str(w_s), ']']})
ylabel('$\bar{R}_t$')
xlabel('Time, $t$ (days)')

figure(2)
plot(days, I)
title('Simulated incidence data')
ylabel('Incidence')
xlabel('Time, $t$ (days)')

figure(3)
plot(0:7, w_s)

