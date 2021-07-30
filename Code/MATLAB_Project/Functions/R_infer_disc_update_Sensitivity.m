function [I, Mean] = R_infer_disc_update_Sensitivity(Serial_Estimate, I_Generation_Method, Hybrid, para_cell)

%This function evaluates the output variables using Bayesian inference
%techniques. The key difference is how the data is generated and how
%accurately the serial interval is tracked.

%It is important to note that whenever there is a switch at time t, we only
%implement the change from time t+1.

%Un-pack
% seed = para.seed; %1
% total_time = para.total_time; %2
% w_s_all_actual = para.w_s_all_actual; %3 %All serial intervals
% w_s_all_recorded = para.w_s_all_recorded; %4
% switch_behaviour = para.switch_behaviour; %5 %These are the times
% update_behaviour = para.update_behaviour; %6 %These are the times
% tau = para.tau; %7
% a = para.a; %8
% b = para.b; %9
% I_0 = para.I_0; %10

rng(para_cell{1}) %Gives reproducible results

%Pre-allocate. We store as matrices because these are changing through
%time.

w_s_actual = zeros(para_cell{2}+1, size(para_cell{3}, 2)); %Don't technically use serial interval for Day 0
w_s_recorded = zeros(para_cell{2}+1, size(para_cell{4}, 2));


w_s_recorded(1:para_cell{6}(1)+1, :) = repmat(para_cell{4}(1, :), para_cell{6}(1)+1, 1);
w_s_recorded(para_cell{6}(end)+1:end, :) = repmat(para_cell{4}(end, :), para_cell{2}+1-para_cell{6}(end), 1);

w_s_actual(1:para_cell{5}(1)+1, :) = repmat(para_cell{3}(1, :), para_cell{5}(1)+1, 1);
w_s_actual(para_cell{5}(end)+1:end, :) = repmat(para_cell{3}(end, :), para_cell{2}+1-para_cell{5}(end), 1);


%THIS WORKS
% w_s_actual = [repmat(para_cell{3}(1, :), para_cell{5}(1)+1, 1); repmat(para_cell{3}(2, :), para_cell{2}-para_cell{5}(1), 1)];
% w_s_recorded = [repmat(para_cell{3}(1, :), para_cell{5}(1)+1, 1); repmat(para_cell{3}(2, :), para_cell{2}-para_cell{5}(1), 1)];

%If we are using hybrid SIs, then both actual and recorded SIs should be
%altered

if isequal(Hybrid, 'Hybrid') || isequal(Hybrid, 'Hybrid-Generation')
    
    if length(para_cell{5})+length(para_cell{6}) ~= 2
        
        disp('Error! Wrong no. of switches/updates for this model version!')
        
        %Could also check that switch times are OK (i.e. equal)
    
    elseif para_cell{5} ~= para_cell{6}
        
        disp('Error! Switch times should equal update times for this model version!')
    
    else
        
        N_o = find(para_cell{3}(1, :), 1, 'last')-1; %Finds last non zero element but takes off 1 because the days start from 0
        
        N_a = find(para_cell{3}(2, :), 1, 'last')-1;
        
        num_hybrids = max(N_a, N_o - N_a);
        
        w_s_hy = zeros(num_hybrids, N_o + 1);
        
        for i = 1:num_hybrids  %doesnt includes the hybrid that is identical to old
            
            %If N_a is larger than N_o-i, then
            
            length_hy = max(N_o-i+1, N_a+1);
            
            tmp = zeros(1, N_o+1); %Length of hybrid
            
            %     if N_o-i >= N_a+1 %The backtracking hasn't reached the new threshold
            
            tmp(1:length_hy) = para_cell{3}(1, 1:length_hy);
            
            if i <= N_a+1
                
                tmp(1:i) = para_cell{3}(1, 1:i);
                
            elseif i > N_a+1
                
                tmp(1:N_a+1) = para_cell{3}(2, 1:N_a+1);
                
            end
            
            tmp = tmp/(sum(tmp)); %Normalization
            
            w_s_hy(i, :) = tmp;
            
        end
        
        if isequal(Hybrid, 'Hybrid')
        
            w_s_recorded(para_cell{6}+1:para_cell{6}+num_hybrids, :) = w_s_hy;
        
        end
        
        w_s_actual(para_cell{5}+1:para_cell{5}+num_hybrids, :) = w_s_hy;
        
    end
        
    %Begin hybridising after the switch times... Maybe check the previous
    %one was OK
    
end
    
% if isequal(Hybrid, 'Non-Hybrid') && length(para_cell{5})>=2 && length(para_cell{6})>=2 
%     
%     for i = 2:length(para_cell{5})
%         
%         w_s_actual(para_cell{5}(i-1)+1:para_cell{5}(i), :) = repmat(para_cell{3}(i, :),para_cell{5}(i) - para_cell{5}(i-1) , 1);
%         
%     end
%     
%     
%     for i = 2:length(para_cell{6})
%         
%         w_s_recorded(para_cell{6}(i-1)+1:para_cell{6}(i), :) = repmat(para_cell{4}(i, :), para_cell{6}(i) - para_cell{6}(i-1), 1);
%         
%     end
%     
% end

% w_s_actual(:, 1) = []; %Reason why we delete is because for the generation, we want to include the most recent day in the genereation of incidence for the present day.
%In contrast, when we infer we want to ignore the present day in inferring
%how the data was generated.

%No need to delete if the first element is not 0

if isequal(I_Generation_Method, 'Trivial')

    R_t = para_celll{end};
    
    I = para_cell{10};
    
    for t = 1:para_cell{2}
        
        I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s_actual(t+1, :)));
        
        I = [I, I_new];
        
    end
        
elseif isequal(I_Generation_Method, 'Variable')
    
    R_t = para_cell{end}([1:1:para_cell{2}+1]);

%     w_s_o(1) = [];
    
    I = para_cell{10};
    
    for t = 1:para_cell{2}
        
        I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_actual(t+1, :)));
        
        I = [I, I_new];
        
    end
    
elseif isequal(I_Generation_Method, 'Data')
    
    I = para_extra.Incidence;
    
end


if isequal(Serial_Estimate, 'Perfect') %Remember that Perfect isnt really perfect here
 
    Shape = zeros(1, para_cell{2});
    Scale = zeros(1, para_cell{2});
    Mean = zeros(1, para_cell{2});
    
    for t = para_cell{7}+1:para_cell{2} % We start from para_cell{7}+1 because the EpiEstim app does too (starting from using the second data point since we cant infer alot from the first)
        
        Shape(t) = para_cell{8} + sum(I(t-para_cell{7}+1:t));
        
        %Calculate summation of Lambdas
        
        for k = t:-1:t-para_cell{7}+1 % k  = t-para_cell{7}+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(k , :)]);
            
        end
        
        I(1) = para_cell{10}; %First case is imported
        
        Scale(t) = 1/(Scale(t)+(1/para_cell{9}));
        
        Mean(t) = Scale(t)*Shape(t);
        
    end
       
elseif isequal(Serial_Estimate, 'Fixed')
    
    Shape = zeros(1, para_cell{2});
    Scale = zeros(1, para_cell{2});
    Mean = zeros(1, para_cell{2});
    
    for t = para_cell{7}+1:para_cell{2}+1
        
        Shape(t) = para_cell{8} + sum(I(t-para_cell{7}+1:t));
        Scale(t) = 0;
        
        %Calculate summation of Lambdas
        
        for k  = t-para_cell{7}+1:t
            
            I_relevant = I(1:k);
            
            Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(1, :)]); %WHY IS THIS ONE DIFFERENT TO SCALE ON PERFECT?
            
        end
        
        I(1) = para_cell{10}; %First case is imported
        
        Scale(t) = 1/(Scale(t)+(1/para_cell{9}));
        
        Mean(t) = Scale(t)*Shape(t);
        
    end

elseif isequal(Serial_Estimate, 'Gradual')

end

end