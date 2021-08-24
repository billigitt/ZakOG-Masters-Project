function [R] = r_to_R(r,w)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(w);

summation_vec = zeros(1, N);

for i = 1:N
   
    summation_vec(i) = w(i)*(exp(-r*(i-0.5)) - exp(-r*(i+0.5)))
    
end

R = r/sum(summation_vec);

end

