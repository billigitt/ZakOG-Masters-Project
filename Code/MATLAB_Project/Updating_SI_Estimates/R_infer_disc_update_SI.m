function [w_s_actual, w_s_recorded, I, Shape, Scale, Mean, Upper, Lower] = R_infer_disc_update_SI(Serial_Estimate, I_Generation_Method, para, para_extra)

%This function evaluates the output variables using Bayesian inference
%techniques. The key difference is how the data is generated and how
%accurately the serial interval is tracked.

%Un-pack
seed = para.seed;
total_time = para.total_time;
w_s_all_actual = para.w_s_all_actual; %All serial intervals
w_s_all_recorded = para.w_s_all_recorded;
switch_behaviour = para.switch_behaviour;
update_behaviour = para.update_behaviour;
tau = para.tau;
a = para.a;
b = para.b;
I_0 = para.I_0;

rng(seed) %Gives reproducible results

w_s_actual = zeros(total_time+1, size(w_s_all_actual, 2));
w_s_recorded = zeros(total_time+1, size(w_s_all_recorded, 2));

disp(size(w_s_recorded(update_behaviour(end)+1:end)))
disp(size(repmat(w_s_all_recorded(end, :), 1, total_time+2-update_behaviour(end))))

w_s_recorded(1:update_behaviour(1)-1, :) = repmat(w_s_all_recorded(1, :), update_behaviour(1)-1, 1);
w_s_recorded(update_behaviour(end):end, :) = repmat(w_s_all_recorded(end, :), total_time+2-update_behaviour(end), 1);

w_s_actual(1:switch_behaviour(1)-1, :) = repmat(w_s_all_actual(1, :), switch_behaviour(1)-1, 1);
w_s_actual(switch_behaviour(end):end, :) = repmat(w_s_all_actual(end, :), total_time+2-switch_behaviour(end), 1);

for i = 2:length(switch_behaviour)

    w_s_actual(switch_behaviour(i-1):switch_behaviour(i)-1, :) = repmat(w_s_all_actual(i, :),switch_behaviour(i) - switch_behaviour(i-1) , 1);
    
end

w_s_actual(:, 1) = [];

for i = 2:length(update_behaviour)

    w_s_recorded(update_behaviour(i-1):update_behaviour(i)-1, :) = repmat(w_s_all_recorded(i, :), update_behaviour(i) - update_behaviour(i-1), 1);
    
end

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
    
end

%%

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