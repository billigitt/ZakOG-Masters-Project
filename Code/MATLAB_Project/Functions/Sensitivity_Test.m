clc
clear all
close all

%% Start with analysis on how the serial inteval affects R_t inference

Key = {'w_s_all_actual', 'w_s_all_recorded'};

N_o = 15;

mu_vec = 2:0.1:13;

% pd_o = zeros(length(mu_vec));
x = linspace(1,N_o,N_o);

w_s_o = zeros(length(mu_vec), N_o);

for i = 1:length(mu_vec)

    pd_o(i) = makedist('Normal','mu',mu_vec(i),'sigma',2);
    
    trunc_o(i) = truncate(pd_o(i), 0, inf);

    w_s_o(i, :) = pdf(trunc_o(i),x)/sum(pdf(trunc_o(i),x));
        
end

pd_a1 = pd_o;

trunc_a1 = trunc_o;

w_s_a1 = w_s_o;

%Get 'vector' of w_s_o and w_s_a1 and encode this into the key parameter.
%Then try to run this.

%Look at error and try to find way of resolving- it is currently trying to
%equate a sclalar with vector... may need to trace which elemtn of
%sensitivity we are on..!

% Epidemiological and Inference parameters

total_time = 100;

tau = 7;

switch_behaviour = 40;

delay = 0;

update_behaviour = switch_behaviour + delay;

I_0 = 100;

para_o = struct('seed', 1, 'total_time', total_time, 'w_s_all_actual', w_s_o, 'w_s_all_recorded', w_s_a1, 'switch_behaviour', switch_behaviour, 'update_behaviour', update_behaviour, 'tau', tau, 'a', 1, 'b', 5, 'I_0', I_0);

R_start = 10;

R_end = 0.6;

days = 0:1:total_time(1);

para_Linear_Vary = struct('R_t',R_start + (R_end(1)-R_start(1))*days/total_time(1));

[Mean_Dif, Area_Dif, para_new] = Sensitivity_Analysis(Key, para_o, para_Linear_Vary, 'Perfect', 'Variable', 'Non-Hybrid')

imagesc(mu_vec, mu_vec, Area_Dif)

colorbar