clc
clear all
close all

set(0,'DefaultFigureColor',[1 1 1])
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaultlinelinewidth', 2)
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex')

%%
%TrivialEstimate_2

R_t = 2; %Fixed R_t

w_s = [0.1 0.2 0.3 0.2 0.1 0.05 0.03 0.02]; %Serial interval, e.g. odds
...of infecting after 2 days is 0.2.

% w_s = [1 0 0 0 0 0];

total_time = 40; %Therefore total time+1 total pieces of data, since I_0
...is at time = 0.

days = 0:total_time;

tau = 5; %time that we sample over to get R_t estimate

%Gamma distribution parameters

a = 1; %This is by solving for mean=5 and stdev=5
b = 5;

%Initial incidence
I_0 = 1;
I = I_0;

%Generate incidence data

%Incidence_Generator_2 takes the whole time series data for I, extracts the
%most recent (last 8 days in this case) incidence reports and then computes
%the dot product with the serial. [The serial is reversed before taking the
%dot product]

for t = 1:total_time
   
    I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s));
%     I_new = R_t*I(end); %Try this to show that current method does not
%     work as gives wrong R_t inference with deterministic data!
    I = [I, I_new];
    
end

Shape = zeros(1, total_time+1);
Scale = zeros(1, total_time+1);
Mean = zeros(1, total_time+1);
Upper = zeros(1, total_time+1);
Lower = zeros(1, total_time+1);

for t = tau+1:total_time+1 %tau+1 (Day tau) is the first you can start 
    ...from. t is the index, NOT the day.
   
    Shape(t) = a + sum(I(t-tau:t)); %shape param at time t = i
    Scale(t) = 0;
    
    for k  = t-tau:t %Index t-tau to t, i.e. day t-tau-1 to t-1
        
        I_switch = I(1:k); %definitely want from k down, may want all
        ...the way down to 1
        
        Scale(t) = Scale(t) + Incidence_Generator_2(I_switch, w_s); %Takes
        ...I_k, I_{k-1} down to length of Serial

    end
    
    Scale(t) = 1/(Scale(t)+(1/b));
    
    Mean(t) = Scale(t)*Shape(t); %Mean of gamma distn is product of two 
    ...pars
    
    Upper(t) = gaminv(0.995, Shape(t), Scale(t));
    
    Lower(t) = gaminv(0.005, Shape(t), Scale(t));
    
end

figure(1)
% plot(days, upper)
% hold on
% plot(days, lower)
plot(days, Mean)

hold on
daysflip = [days, fliplr(days)];
inBetween = [Lower, fliplr(Upper)];
I_switch = fill(daysflip, inBetween, 'r', 'FaceAlpha', 0.25, ...
    'edgealpha', 0);
title({['Trivial $R_t$ (=', num2str(R_t), ') inference with'];...
    ['$w_s = [$',num2str(w_s), ']']})
ylabel('$\bar{R}_t$')
xlabel('Time, $t$ (days)')

figure(2)
plot(days, I)
hold on
plot(days, R_t.^days, 'r')
title('Simulated incidence data')
ylabel('Incidence')
xlabel('Time, $t$ (days)')

% figure(3)
% plot(0:7, w_s, 'k')
% hold on
% plot(days, R_t.^days, 'r')
