function [w_s_actual, w_s_recorded, I, Shape, Scale, Mean, Upper, Lower] = R_infer_disc_multiupdate_SI(Serial_Estimate, I_Generation_Method, Hybrid, para, para_extra)

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

w_s_actual = zeros(total_time+1, size(w_s_all_actual, 2));
w_s_recorded = zeros(total_time+1, size(w_s_all_recorded, 2));

disp(size(w_s_recorded(update_behaviour(end)+1:end)))
disp(size(repmat(w_s_all_recorded(end, :), 1, total_time+2-update_behaviour(end))))

w_s_recorded(1:update_behaviour(1), :) = repmat(w_s_all_recorded(1, :), update_behaviour(1), 1);
w_s_recorded(update_behaviour(end)+1:end, :) = repmat(w_s_all_recorded(end, :), total_time+1-update_behaviour(end), 1);

w_s_actual(1:switch_behaviour(1), :) = repmat(w_s_all_actual(1, :), switch_behaviour(1), 1);
w_s_actual(switch_behaviour(end)+1:end, :) = repmat(w_s_all_actual(end, :), total_time+1-switch_behaviour(end), 1);

%If we are using hybrid SIs, then both actual and recorded SIs should be
%altered

if isequal(Hybrid, 'Hybrid') || isequal(Hybrid, 'Hybrid-Generation')
    
    if switch_behaviour ~= update_behaviour
        
        disp('Error! Switch times should equal update times for this model version!')
        
    else
        
        %Generate vector of size of serial intyervals, N
        
        num_serials = size(w_s_all_actual, 1);
        
        N = zeros(num_serials, 1);
        
        for r = 1:num_serials %No. of serial intervals
           
            N(r) = find(w_s_all_actual(r, :), 1, 'last')-1;
            
        end
        
        %% Now we need to do the hybridisation for each transition point.
        
        for w = 1:length(switch_behaviour)
            
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
                
                w_s_hy(i, :, w) = tmp;
                
            end
            
        end
        
        if isequal(Hybrid, 'Hybrid')
            
            w_s_recorded(switch_behaviour+1:switch_behaviour+num_hybrids, :) = w_s_hy;
            
        end
        
        w_s_actual(switch_behaviour+1:switch_behaviour+num_hybrids, :) = w_s_hy;
        
    end
        
    %Begin hybridising after the switch times... Maybe check the previous
    %one was OK
    
end
    
if isequal(Hybrid, 'Non-Hybrid') && length(switch_behaviour)>=2 && length(update_behaviour)>=2 
    
    for i = 2:length(switch_behaviour)
        
        w_s_actual(switch_behaviour(i-1):switch_behaviour(i)-1, :) = repmat(w_s_all_actual(i, :),switch_behaviour(i) - switch_behaviour(i-1) , 1);
        
    end
    
    
    for i = 2:length(update_behaviour)
        
        w_s_recorded(update_behaviour(i-1):update_behaviour(i)-1, :) = repmat(w_s_all_recorded(i, :), update_behaviour(i) - update_behaviour(i-1), 1);
        
    end
    
end

w_s_actual(:, 1) = [];

if isequal(I_Generation_Method, 'Trivial')

    R_t = para_extra.R_t;
    
    I = I_0;
    
    for t = 1:total_time
        
        I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s_actual(t, :)));
        
        I = [I, I_new];
        
    end
        
elseif isequal(I_Generation_Method, 'Variable')
    
    R_t = para_extra.R_t([0:1:total_time]);

%     w_s_o(1) = [];
    
    I = I_0;
    
    for t = 1:total_time
        
        I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_actual(t, :)));
        
        I = [I, I_new];
        
    end
    
elseif isequal(I_Generation_Method, 'Data')
    
    I = para_extra.Incidence;
    
end


if isequal(Serial_Estimate, 'Perfect') %Remember that Perfect isnt really perfect here
 
    Shape = zeros(1, total_time+1);
    Scale = zeros(1, total_time+1);
    Mean = zeros(1, total_time+1);
    Upper = zeros(1, total_time+1);
    Lower = zeros(1, total_time+1);
    
    for t = tau+1:total_time+1
        
        Shape(t) = a + sum(I(t-tau+1:t));
        Scale(t) = 0;
        
        %Calculate summation of Lambdas
        
        for k  = t-tau+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [w_s_recorded(t, :)]);
            
        end
        
        I(1) = I_0; %First case is imported
        
        Scale(t) = 1/(Scale(t)+(1/b));
        
        Mean(t) = Scale(t)*Shape(t);
        
        Upper(t) = gaminv(0.975, Shape(t), Scale(t));
        
        Lower(t) = gaminv(0.025, Shape(t), Scale(t));
        
    end
       
elseif isequal(Serial_Estimate, 'Fixed')
    
    Shape = zeros(1, total_time+1);
    Scale = zeros(1, total_time+1);
    Mean = zeros(1, total_time+1);
    Upper = zeros(1, total_time+1);
    Lower = zeros(1, total_time+1);
    
    
    for t = tau+1:total_time+1
        
        Shape(t) = a + sum(I(t-tau+1:t));
        Scale(t) = 0;
        
        %Calculate summation of Lambdas
        
        for k  = t-tau+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(1, :)]);
            
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

