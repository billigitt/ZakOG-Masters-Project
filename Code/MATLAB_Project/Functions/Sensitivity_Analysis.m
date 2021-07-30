function [Mean_Dif, Area_Dif, para_all] = Sensitivity_Analysis(key_para_string, para, para_extra, Serial_Estimate, I_Generation_Method, Hybrid)

%This function performs sensitivity analyses on two key parameters and
%determines both the end time offset of the unupdated inference versus the
%updated inference, as well as the cumulative offset (area between incidences)

%This should illustrate how these two key parameters interact.

%para is all of the parameters, and two of the fields within para will be
%vectors (or matrices if the serial intervals are chosen)

%Plan

%Find which one is the key parameter

%On each for loop change this value
%Massive if loop for if key_string = 'total_time' else if key_string =
%'switch' etc

%Un-pack

 %Gives reproducible results

%Identify key parameter

para_cell = struct2cell(para);
para_extra_cell = struct2cell(para_extra);


para_all = [para_cell; para_extra_cell];
FieldNames_1 = fieldnames(para);
FieldNames_2 = fieldnames(para_extra);

FieldNames = [FieldNames_1; FieldNames_2];

idx1 = find(contains(FieldNames,key_para_string(1)));
idx2 = find(contains(FieldNames,key_para_string(2)));

length_1 = size(para_cell{idx1}, 2); %2 because this will work generally, including for the serial intervals

if idx2 > length(FieldNames_1)

    length_2 = size(para_extra_cell{1}, 2);

else 
    
    length_2 = size(para_cell{idx2}, 2);
    
end

% seed = para.seed;
% total_time = para.total_time;
% w_s_all_actual = para.w_s_all_actual; %All serial intervals
% w_s_all_recorded = para.w_s_all_recorded;
% switch_behaviour = para.switch_behaviour; %These are the times
% update_behaviour = para.update_behaviour; %These are the times
% tau = para.tau;
% a = para.a;
% b = para.b;
% I_0 = para.I_0;

%Get string that is asspciated to index
%Find way to get variable of string
% NEED TO UPDATE R_T FOR DIFFERENT LENGTHS OF TIME

para_tmp = para_all;

Mean_Dif = zeros(length_1, length_2);

Area_Dif = zeros(length_1, length_2);

rng(para_tmp{1})

for index_1 = 1:length_1
    
    para_tmp{idx1} = para_all{idx1}(index_1);
    
    %Update key parameter_1
    
    for index_2 = 1:length_2
       
        para_tmp{idx2} = para_all{idx2}(index_2);
        
        [~, tmp_Mean_Update] = R_infer_disc_update_Sensitivity(Serial_Estimate, I_Generation_Method, Hybrid, para_tmp);
        
        para_tmp{4} = para_tmp{3};
        
        [~, tmp_Mean_Non_Update] = R_infer_disc_update_Sensitivity(Serial_Estimate, I_Generation_Method, Hybrid, para_tmp);
        
        mean_NU = tmp_Mean_Non_Update(para_tmp{5}+1:end);
        
        mean_U = tmp_Mean_Update(para_tmp{5}+1:end);
        
        Mean_Dif(index_1, index_2) = mean_NU(end)-mean_U(end);
        
        Area_Dif(index_1, index_2) = trapz(para_tmp{5}+1:para_tmp{2}, abs(mean_NU-mean_U));
        
    end
    
end

%Get length of these simulations

% length_1 = 
% 
% length_2 = 

% if idx2 > length(FieldNames_1) %R_t is one of the key parameters
%    
%     String_2 = 'R_t';
%     
%     String_1 = FieldNames(idx1);
%     
% else
%     
%     String_1 = Field
%     
% end




%% Comment from here

% %Pre-allocate. We store as matrices because these are changing through
% %time.
% 
% w_s_actual = zeros(total_time+1, size(w_s_all_actual, 2)); %Don't technically use serial interval for Day 0
% w_s_recorded = zeros(total_time+1, size(w_s_all_recorded, 2));
% 
% 
% w_s_recorded(1:update_behaviour(1)+1, :) = repmat(w_s_all_recorded(1, :), update_behaviour(1)+1, 1);
% w_s_recorded(update_behaviour(end)+1:end, :) = repmat(w_s_all_recorded(end, :), total_time+1-update_behaviour(end), 1);
% 
% w_s_actual(1:switch_behaviour(1)+1, :) = repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1);
% w_s_actual(switch_behaviour(end)+1:end, :) = repmat(w_s_all_actual(end, :), total_time+1-switch_behaviour(end), 1);
% 
% 
% %THIS WORKS
% % w_s_actual = [repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1); repmat(w_s_all_actual(2, :), total_time-switch_behaviour(1), 1)];
% % w_s_recorded = [repmat(w_s_all_actual(1, :), switch_behaviour(1)+1, 1); repmat(w_s_all_actual(2, :), total_time-switch_behaviour(1), 1)];
% 
% %If we are using hybrid SIs, then both actual and recorded SIs should be
% %altered
% 
% 
% if isequal(I_Generation_Method, 'Trivial')
% 
%     R_t = para_extra.R_t;
%     
%     I = I_0;
%     
%     for t = 1:total_time
%         
%         I_new = poissrnd(R_t*Incidence_Generator_2(I, w_s_actual(t+1, :)));
%         
%         I = [I, I_new];
%         
%     end
%         
% elseif isequal(I_Generation_Method, 'Variable')
%     
%     R_t = para_extra.R_t([1:1:total_time+1]);
% 
% %     w_s_o(1) = [];
%     
%     I = I_0;
%     
%     for t = 1:total_time
%         
%         I_new = poissrnd(R_t(t)*Incidence_Generator_2(I, w_s_actual(t+1, :)));
%         
%         I = [I, I_new];
%         
%     end
%     
% elseif isequal(I_Generation_Method, 'Data')
%     
%     I = para_extra.Incidence;
%     
% end
% 
% 
% if isequal(Serial_Estimate, 'Perfect') %Remember that Perfect isnt really perfect here
%  
%     Shape = zeros(1, total_time);
%     Scale = zeros(1, total_time);
%     Mean = zeros(1, total_time);
%     Upper = zeros(1, total_time);
%     Lower = zeros(1, total_time);
%     
%     for t = tau+1:total_time % We start from tau+1 because the EpiEstim app does too (starting from using the second data point since we cant infer alot from the first)
%         
%         Shape(t) = a + sum(I(t-tau+1:t));
%         
%         %Calculate summation of Lambdas
%         
%         for k = t:-1:t-tau+1 % k  = t-tau+1:t
%             
%             I_relevant = I(1:k);
%             
%             Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(k , :)]);
%             
%         end
%         
%         I(1) = I_0; %First case is imported
%         
%         Scale(t) = 1/(Scale(t)+(1/b));
%         
%         Mean(t) = Scale(t)*Shape(t);
%         
%         Upper(t) = gaminv(0.975, Shape(t), Scale(t));
%         
%         Lower(t) = gaminv(0.025, Shape(t), Scale(t));
%         
%     end
%        
% elseif isequal(Serial_Estimate, 'Fixed')
%     
%     Shape = zeros(1, total_time);
%     Scale = zeros(1, total_time);
%     Mean = zeros(1, total_time);
%     Upper = zeros(1, total_time);
%     Lower = zeros(1, total_time);
%     
%     
%     for t = tau+1:total_time+1
%         
%         Shape(t) = a + sum(I(t-tau+1:t));
%         Scale(t) = 0;
%         
%         %Calculate summation of Lambdas
%         
%         for k  = t-tau+1:t
%             
%             I_relevant = I(1:k);
%             
%             Scale(t) = Scale(t) + Incidence_Generator_2(I_relevant, [0 w_s_recorded(1, :)]); %WHY IS THIS ONE DIFFERENT TO SCALE ON PERFECT?
%             
%         end
%         
%         I(1) = I_0; %First case is imported
%         
%         Scale(t) = 1/(Scale(t)+(1/b));
%         
%         Mean(t) = Scale(t)*Shape(t);
%         
%         Upper(t) = gaminv(0.975, Shape(t), Scale(t));
%         
%         Lower(t) = gaminv(0.025, Shape(t), Scale(t));
%         
%     end
% 
% elseif isequal(Serial_Estimate, 'Gradual')
% 
% end
% 
% end