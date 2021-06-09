function Lambda = Incidence_Generator(I_data, Serial) %
    
    length_I = length(I_data);
    
    length_S = length(Serial);

    if length_I< length_S
       
        I_recent = [I_data, zeros(1,length_S - length_I)];
        
    else
        
        I_recent = I_data(1:length_S);
        
    end

    Lambda = dot(I_recent, Serial);

end