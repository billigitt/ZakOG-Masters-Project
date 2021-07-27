function [w_s_actual, w_s_recorded, I, Shape, Scale, Mean, Upper, Lower] = R_infer_disc_update_SI(Serial_Estimate, I_Generation_Method, Hybrid, para, para_extra)

%This function evaluates the output variables using Bayesian inference
%techniques. The key difference is how the data is generated and how
%accurately the serial interval is tracked.

%It is important to note that whenever there is a switch at time t, we only
%implement the change from time t+1.

%Un-pack
seed = para.seed;
total_time = para.total_time;
w_s_all_actual = para.w_s_all_actual; %All serial intervals
w_s_all_recorded = para.w_s_all_recorded;
switch_behaviour = para.switch_behaviour; %These are the times
update_behaviour = para.update_behaviour; %These are the times
tau = para.tau;
a = para.a;
b = para.b;
I_0 = para.I_0;

rng(seed) %Gives reproducible results

%Pre-allocate. We store as matrices because these are changing through
%time.

w_s_actual = zeros(total_time+1, size(w_s_all_actual, 2)); %Don't technically use serial interval for Day 0
w_s_recorded = zeros(total_time+1, size(w_s_all_recorded, 2));


w_s_recorded(1:update_behaviour(1)+1, :) = repmat(w_s_all_recorded(1, :), update_behaviour(1)+1, 1);
w_s_recorded(update_behaviour(end)+1:end, :) = repmat(w_s_all_recorded(end, :), total_time+1-update_behaviour(end), 1);

w_s_actual(1:switch_behaviour(1)+1, :) = repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1);
w_s_actual(switch_behaviour(end)+1:end, :) = repmat(w_s_all_actual(end, :), total_time+1-switch_behaviour(end), 1);


%THIS WORKS
% w_s_actual = [repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1); repmat(w_s_all_actual(2, :), total_time-switch_behaviour(1), 1)];
% w_s_recorded = [repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1); repmat(w_s_all_actual(2, :), total_time-switch_behaviour(1), 1)];

%If we are using hybrid SIs, then both actual and recorded SIs should be
%altered

if isequal(Hybrid, 'Hybrid') || isequal(Hybrid, 'Hybrid-Generation')
    
    if length(switch_behaviour)+length(update_behaviour) ~= 2
        
        disp('Error! Wrong no. of switches/updates for this model version!')
        
        %Could also check that switch times are OK (i.e. equal)
    
    elseif switch_behaviour ~= update_behaviour
        
        disp('Error! Switch times should equal update times for this model version!')
    
    else
        
        N_o = find(w_s_all_actual(1, :), 1, 'last')-1; %Finds last non zero element but takes off 1 because the days start from 0
        
        N_a = find(w_s_all_actual(2, :), 1, 'last')-1;
        
        num_hybrids = max(N_a, N_o - N_a);
        
        w_s_hy = zeros(num_hybrids, N_o + 1);
        
        for i = 1:num_hybrids  %doesnt includes the hybrid that is identical to old
            
            %If N_a is larger than N_o-i, then
            
            length_hy = max(N_o-i+1, N_a+1);
            
            tmp = zeros(1, N_o+1); %Length of hybrid
            
            %     if N_o-i >= N_a+1 %The backtracking hasn't reached the new threshold
            
            tmp(1:length_hy) = w_s_all_actual(1, 1:length_hy);
            
            if i <= N_a+1
                
                tmp(1:i) = w_s_all_actual(1, 1:i);
                
            elseif i > N_a+1
                
                tmp(1:N_a+1) = w_s_all_actual(2, 1:N_a+1);
                
            end
            
            tmp = tmp/(sum(tmp)); %Normalization
            
            w_s_hy(i, :) = tmp;
            
        end
        
        if isequal(Hybrid, 'Hybrid')
        
            w_s_recorded(update_behaviour+1:update_behaviour+num_hybrids, :) = w_s_hy;
        
        end
        
        w_s_actual(switch_behaviour+1:switch_behaviour+num_hybrids, :) = w_s_hy;
        
    end
        
    %Begin hybridising after the switch times... Maybe check the previous
    %one was OK
    
end
    
% if isequal(Hybrid, 'Non-Hybrid') && length(switch_behaviour)>=2 && length(update_behaviour)>=2 
%     
%     for i = 2:length(switch_behaviour)
%         
%         w_s_actual(switch_behaviour(i-1)+1:switch_behaviour(i), :) = repmat(w_s_all_actual(i, :),switch_behaviour(i) - switch_behaviour(i-1) , 1);
%         
%     end
%     
%     
%     for i = 2:length(update_behaviour)
%         
%         w_s_recorded(update_behaviour(i-1)+1:update_behaviour(i), :) = repmat(w_s_all_recorded(i, :), update_behaviour(i) - update_behaviour(i-1), 1);
%         
%     end
%     
% end

% w_s_actual(:, 1) = []; %Reason why we delete is because for the generation, we want to include the most recent day in the genereation of incidence for the present day.
%In contrast, when we infer we want to ignore the present day in inferring
%how the data was generated.

%No need to delete if the first element is not 0

if isequal(I_Generation_Method, 'Trivial')

    R_t = para_extra.R_t;
    
    I = I_0;
    
    for t = 1:total_time
        
        I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s_actual(t+1, :)));
        
        I = [I, I_new];
        
    end
        
elseif isequal(I_Generation_Method, 'Variable')
    
    R_t = para_extra.R_t([1:1:total_time+1]);

%     w_s_o(1) = [];
    
    I = I_0;
    
    for t = 1:total_time
        
        I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_actual(t+1, :)));
        
        I = [I, I_new];
        
    end
    
elseif isequal(I_Generation_Method, 'Data')
    
    I = para_extra.Incidence;
    
end


if isequal(Serial_Estimate, 'Perfect') %Remember that Perfect isnt really perfect here
 
    Shape = zeros(1, total_time);
    Scale = zeros(1, total_time);
    Mean = zeros(1, total_time);
    Upper = zeros(1, total_time);
    Lower = zeros(1, total_time);
    
    for t = tau+1:total_time % We start from tau+1 because the EpiEstim app does too (starting from using the second data point since we cant infer alot from the first)
        
        Shape(t) = a + sum(I(t-tau+1:t));
        
        %Calculate summation of Lambdas
        
        for k = t:-1:t-tau+1 % k  = t-tau+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(k , :)]);
            
        end
        
        I(1) = I_0; %First case is imported
        
        Scale(t) = 1/(Scale(t)+(1/b));
        
        Mean(t) = Scale(t)*Shape(t);
        
        Upper(t) = gaminv(0.975, Shape(t), Scale(t));
        
        Lower(t) = gaminv(0.025, Shape(t), Scale(t));
        
    end
       
elseif isequal(Serial_Estimate, 'Fixed')
    
    Shape = zeros(1, total_time);
    Scale = zeros(1, total_time);
    Mean = zeros(1, total_time);
    Upper = zeros(1, total_time);
    Lower = zeros(1, total_time);
    
    
    for t = tau+1:total_time+1
        
        Shape(t) = a + sum(I(t-tau+1:t));
        Scale(t) = 0;
        
        %Calculate summation of Lambdas
        
        for k  = t-tau+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(1, :)]); %WHY IS THIS ONE DIFFERENT TO SCALE ON PERFECT?
            
        end
        
        I(1) = I_0; %First case is imported
        
        Scale(t) = 1/(Scale(t)+(1/b));
        
        Mean(t) = Scale(t)*Shape(t);
        
        Upper(t) = gaminv(0.975, Shape(t), Scale(t));
        
        Lower(t) = gaminv(0.025, Shape(t), Scale(t));
        
    end

elseif isequal(Serial_Estimate, 'Gradual')

end

end