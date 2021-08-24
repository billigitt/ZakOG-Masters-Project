function [delta, w_a] = Delta_Generator(w_o, mu_d, var_d)

N = length(w_o);

delta = zeros(1, N);

A = [N N*(N+1)/2 N*(N+1)*(2*N+1)/6; N*(N+1)/2 N*(N+1)*(2*N+1)/6 (N*(N+1)/2)^2; N*(N+1)*(2*N+1)/6 (N*(N+1)/2)^2 N*(6*N^4 +15*N^3 +10*N^2-1)/30];

u = [0; mu_d; var_d];

k = A\u;

for i = 1:N
   
    delta(i) = dot(k, [1 i, i^2]);
    
end

w_a = delta + w_o;

end