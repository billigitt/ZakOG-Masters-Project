function [w_s_t, I, Shape, Scale, Mean, Upper, Lower] = R_infer_update_SI_2(Serial_Estimate, I_Generation_Method, para, para_extra)

%This function evaluates the output variables using Bayesian inference
%techniques. The key difference is how the data is generated and how
%accurately the serial interval is tracked.

%Un-pack
seed = para.seed;
total_time = para.total_time;
w_s_o = para.w_s_o; %Normal generated beforehand
w_s_f = para.w_s_f; %Truncated normal generated beforehand
tau = para.tau;
a = para.a;
b = para.b;
I_0 = para.I_0;

rng(seed) %Gives reproducible results

w_s_t = zeros(total_time+1, length(w_s_o));

w_s_t(1, :) = w_s_o;


if isequal(I_Generation_Method, 'Trivial')

    R_t = para_extra.R_t;

    w_s_o(1) = [];
    
    I = I_0;
    
    for t = 1:total_time
        
        w_s_t(t+1, :) = (t*w_s_f +(total_time-t)*[0 w_s_o])/total_time;
        
        I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s_t(t, :)));
        
        I = [I, I_new];
        
    end
        
elseif isequal(I_Generation_Method, 'Variable')
    
    R_t = para_extra.R_t([0:1:total_time]);

    w_s_o(1) = [];
    
    I = I_0;
    
    for t = 1:total_time
        
        w_s_t(t+1, :) = (t*w_s_f +(total_time-t)*[0 w_s_o])/total_time;
        
        I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_t(t, :)));
        
        I = [I, I_new];
        
    end
    
elseif isequal(I_Generation_Method, 'Data')
    
end

if isequal(Serial_Estimate, 'Perfect')
 
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
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_t(k, :)]);
            
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
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_o]);
            
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