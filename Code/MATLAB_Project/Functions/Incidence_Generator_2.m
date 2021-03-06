function Gamma = Incidence_Generator_2(I_data, Serial) %
  
%Matches with arranging incidence from I_0 to I_end

    length_I = length(I_data);
    
    length_S = length(Serial);

    %Get the serial and the incidence to be the same length so we can dot
    %product
    
    if length_I < length_S
       
        I_recent = [zeros(1,length_S - length_I), I_data];
        
    else
        
        I_recent = I_data(end-length_S+1:end);
        
    end

    Gamma = dot(I_recent, fliplr(Serial)); %We now have to flip the Serial

end