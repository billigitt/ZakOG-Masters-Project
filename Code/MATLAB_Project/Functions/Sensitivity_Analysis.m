function [Mean_Dif, Area_Dif, para_all] = Sensitivity_Analysis(key_para_string, para, para_extra, Serial_Estimate, I_Generation_Method, Hybrid)

%This function performs sensitivity analyses on two key parameters and
%determines both the end time offset of the unupdated inference versus the
%updated inference, as well as the cumulative offset (area between incidences)

%This should illustrate how these two key parameters interact.

%para is all of the parameters, and two of the fields within para will be
%vectors (or matrices if the serial intervals are chosen)

%Plan

%Find which one is the key parameter (we cannot have idx=4 as a kaye parameter. we just have idx =3 as one)

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

size_1 = size(para_cell{idx1}); % Used to be 2 because this will work generally, including for the serial intervals

length_1 = size_1(end); %Works for 3d matrix too

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

%Create cells of all paremters for updated and non-updated serials

para_tmp = para_all;

para_tmp_NU = para_all;

if idx1 ~= 3

para_tmp_NU{4}(2, :) = para_all{4}(1, :); %Recorded for updated equals the actual 

end

Mean_Dif = zeros(length_1, length_2);

Area_Dif = zeros(length_1, length_2);

rng(para_tmp{1})

idx1

for index_1 = 1:length_1
    
    
    if idx1 == 3 || idx1 == 4 %In this case a 3d matrix
        
        para_tmp{idx1} = para_all{idx1}(:, :, index_1); %Since it will be a 3d matrix

        para_tmp_NU{idx1} = para_all{idx1}(:, :, index_1);
    
    else
    
        para_tmp{idx1} = para_all{idx1}(index_1);

        para_tmp_NU{idx1} = para_all{idx1}(index_1);
        
    end

    
    if idx1 == 3 %In the case when idx1==3, we need to make sure that the recorded is updated for tmp but not for tmp_NU!
        
        para_tmp{4} = para_all{3}(:, :, index_1);
        
        para_tmp_NU{4} = para_all{3}(:, :, index_1);
        
        para_tmp_NU{4}(2, :) = para_all{3}(1, :, index_1); %No update!
        
    end

    %Update key parameter_1
    
    for index_2 = 1:length_2
        
       
        if idx2 == 3 || idx2 == 4 %In this case para_all{idx} is a vector
            
            para_tmp{idx2} = para_all{idx2}(:, :, index_2);
            
            para_tmp_NU{idx2} = para_all{idx2}(:, :, index_2);
            
        else
        
            para_tmp{idx2} = para_all{idx2}(index_2);
            
            para_tmp_NU{idx2} = para_all{idx2}(index_2);
            
        end
        
%         if idx2 == 3
%             
%            para_tmp{4} = para_tmp{3}(:, :, index_2);
%             
%         end
        
        [~, tmp_Mean_Update] = R_infer_disc_update_Sensitivity(Serial_Estimate, I_Generation_Method, Hybrid, para_tmp);
        
        [~, tmp_Mean_Non_Update] = R_infer_disc_update_Sensitivity(Serial_Estimate, I_Generation_Method, Hybrid, para_tmp_NU);
        
        mean_NU = tmp_Mean_Non_Update(para_tmp_NU{5}+1:end);
        
        mean_U = tmp_Mean_Update(para_tmp{5}+1:end);
        
        Mean_Dif(index_1, index_2) = abs((mean_NU(end)-mean_U(end)))/mean_U(end);
        
        Area_Dif(index_1, index_2) = trapz(para_tmp{5}+1:para_tmp{2}, abs(mean_NU-mean_U))/(trapz(para_tmp{5}+1:para_tmp{2}, mean_U));
        
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
