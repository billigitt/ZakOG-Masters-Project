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

%% Empirical evidence for how we can predict the non zero c value for which SI is not important

%Plan:

%Create a vector of mu_d (mean of Delta) running from -1 to 1. From this we
%find generate Delta-vectors. From that, we use fsolve to estimate c-values
%that yield h(c) = 1 where c =/= 0. Then we use the sensitivity analyses
%with key parameters, R and w_s_actual. We input w_s_actual with the
%(delta+w_o) values that we have already found (choosing a suitable w_o)
%and then (once we find a way to get R to correspond to c- should be
%doable?), we can then overlay the c-estimate with the sensitivity analyses
%and we should find they match up perfectly!

%% load('Large_Simulation_1') to skip running code- this will take roughly 8 hours!

Num_Serials = 100;

mu_d = linspace(.3, -.3, Num_Serials); %Differences in means of old and new serials

w_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

N = 10;

Days = 1:N;

Delta_Matrix = zeros(1, N, Num_Serials);

% options = optimset('TolFun',1e-20,'MaxFunEvals',1e8,'Maxiter',1e9);

% C_star = zeros(1, Num_Serials); %Vector of r-values that give h(r) = 1
% 
% R_star = zeros(1, Num_Serials); %Corresponding R-numbers (using eqn 3.6  in Wallinga-Lipsitch)

epsilon = 1e-6; %enables the domain to be roughly (0, x], i.e. non inclusive to 0

w = zeros(1, Num_Serials);

for i = 1:Num_Serials
    
    [Delta_Matrix(1, :, i), ~] = Delta_Generator_Quadratic_1(w_o, mu_d(i), +5); %Generates Deltas with different means. Set differing variance to 5 as enables good plot
    
    w(i) = (w_o(1) + Delta_Matrix(1, 1, i))/w_o(1);
    
%     fun = @(x) h_Transform(x, Delta_Matrix(1, :, i));
    
%     if fun(epsilon)<0 %There is a switch in the sign of mu_d which means the root will be found in a new place. Since fzero relies on a sign change we have to track that appropriatley
%     
%         C_star(i) = fzero(fun,[epsilon 10], options); %If we extend the range of mu_d, we may need larger than 10.
%     
%     else
%     
%         C_star(i) = fzero(fun,[-10 -epsilon], options);
%         
%     end
%         
%     R_star(i) = r_to_R(C_star(i), w_o + Delta_Matrix(1, :, i)); %Converts this c (r epidemic growth value) to an R number, using Wallinga-Lipsitch.
%     
end

%% Sensitivity part

Key = {'w_s_all_actual', 'R_t'};

w_o = [0.08 0.1 0.16 0.14 0.12 0.10 0.09 0.08 0.07 0.06];

w_a_Matrix = repmat(w_o, 1, 1, Num_Serials)+ Delta_Matrix;

w_all_actual = zeros(2, N, Num_Serials);

w_all_actual(1, :, :) = repmat(w_o, 1, 1, Num_Serials);

w_all_actual(2, :, :) = w_a_Matrix;

total_time = 100;

tau = 7;

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

I_0 = 1;

para = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_all_actual, 'w_s_all_recorded', w_all_actual, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);

R_start = 1;

R_end = 6000;
 
days = 0:1:total_time(1);

para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));

Num_Sims = 100;

%%
    
[Ratio, para_new] = Sensitivity_Analysis_H(Key, para, para_Linear_Vary, 'Perfect', 'Trivial_and_Expectation', 'Non-Hybrid');

%% Plots

R_t = para_Linear_Vary.R_t;

%This figure shows the spread of Deltas that we will be studying. They are
%calculated based on what mu_d value we desire.

figure(1)
clf
plot(Days, Delta_Matrix(1, :, 1), 'color', C(1,:))

hold on

for i = 2:20
   
    plot(Days, Delta_Matrix(1, :, i), 'color', C(1, :))
    
end

for i = 21:40
   
    plot(Days, Delta_Matrix(1, :, i), 'color', C(2, :))
    
end

for i = 41:60
   
    plot(Days, Delta_Matrix(1, :, i), 'color', C(3, :))
    
end

for i = 61:80
   
    plot(Days, Delta_Matrix(1, :, i), 'color', C(4, :))
    
end

for i = 81:100
   
    plot(Days, Delta_Matrix(1, :, i), 'color', C(5, :))
    
end



%This figure should tell us when to expect R>1 or R<1  and for what mu
%values. If the output is less than 0, it means that exactly one of 
%Delta_1 or mu_d are less than 0. Otherwise, they both have the same sign.
%The theory tells us that if they have the same sign, the output is
%positive and there will be R values >1 s.t h=1. If the output is negative,
%then the opposite is true and there are c/r values<0 (R<1) s.t h(c)=1.

figure(2)
clf
vec = zeros(1, Num_Serials);

for i = 1:Num_Serials

    vec(i)  = Delta_Matrix(1, 1, i);

end



plot(vec.*mu_d)

%This figure is less useful in the long run compared to R_star

% figure(3)
% clf
% plot(mu_d, C_star)
% 
% figure(4)
% clf
% plot(mu_d, R_star)

figure(5)
clf
imagesc(R_t, mu_d, Ratio)

hold on

% h = plot(R_star, mu_d, 'k', 'LineStyle', '--');
xlabel('True $R_t$')
ylabel('$\mu_{\Delta}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$h(r(R_t))$', 'Interpreter', 'latex')

% legend(h, 'Numerical solve for $v(R(r)) \cdot  \mu_{\Delta}=0$', 'Location', 'SouthEast')
%%

figure(6)
clf
imagesc(R_t, mu_d, Ratio)
colormap parula
hold on

% h = plot(R_star, mu_d, 'k', 'LineStyle', '--');
xlabel('True $R_t$')
ylabel('$\mu_{\mbox{\boldmath $\Delta$}}$')
set(gca,'YDir','normal')

c = colorbar;
set(c,'TickLabelInterpreter','latex')
ylabel(c, '$h(r(R_t))$', 'Interpreter', 'latex')

% legend(h, 'Numerical solve for $0 = \mu_{\Delta} \cdot $\boldmath$v$ $_R$', 'Location', 'SouthEast')

% legend(h, 'Numerical solve for $\mbox{\boldmath $v$}(r(R_t))\cdot\mbox{\boldmath $\Delta$}=0$', 'Location', 'SouthEast')

figure(7)
clf
plot(w, 'color', [0 0 0])
hold on

Num_R = size(Ratio, 2);

error_long_term_estimate = zeros(1, Num_R);

for j = 1:Num_R
   
    error_long_term_estimate(j) = abs(sum(Ratio(:, j)-w')/sum(w'));
    
    plot(Ratio(:, j), 'color', [0 0 0] + j*(C(4, :))/Num_R, 'Linewidth', 0.02)
    
end

figure(8)
clf
X = linspace(0.1, 6000, 100);


m = (log10(0.0440)-log10(4.4*10^-4))/(log10(61)-log10(6000));
c = log10(0.044)+1.0036*log10(61);

h(1) = semilogy(R_t, error_long_term_estimate, 'o', 'color', C(2, :));

hold on


h(2) = semilogy(X, (10^c)*X.^m, 'color', C(4, :));

ylabel('$g(R_t)$')
xlabel('$R_t$')

title('Empirical evidence supporting $\lim_{r \rightarrow \infty} h(r) = \frac{w_1^a}{w_1^o}$')

legend(h, '$g(R_t) = |(\mbox{\boldmath $h$}(R_t)- \mbox{\boldmath $w$})\cdot \mbox{\boldmath $1$}|/|\mbox{\boldmath $w$}|$', '$g(R_t) = 10^{0.435}\cdot x^{-1.00360}$', 'Location', 'NorthWest')

axes('position',[.45 .4 .4 .25])
box on % put box around new pair of axes
indexOfInterest = (X < 400) & (X > 0.01); % range of t near perturbation
semilogy(R_t(indexOfInterest),error_long_term_estimate(indexOfInterest), 'o', 'color', C(2, :)) % plot on new axes
axis tight

hold on
semilogy(X(indexOfInterest), (10^c)*X(indexOfInterest).^m, 'color', C(4, :))




figure(9)
clf
plot(Days, Delta_Matrix(1, :, 1), 'color', C(1,:))


t = linspace(0,2*pi,1000); % values in time
x = cos(t); % a sinusoid
perturbation = 0.1*exp((-(t-5*pi/4).^2)/.01).*sin(200*t); % a perturbation
signal = x+perturbation; % a signal to plot
% plot signal on new figure

figure(10)
clf
plot(t,x+perturbation)
xlabel('time'),ylabel('signal')
% create a new pair of axes inside current figure
axes('position',[.75 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
axis tight