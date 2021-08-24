function [F] = F3D(k, N, mu_d, var_d)

F = zeros(1, 3);

for i = 1:N
   
    F(1) = F(1) + k(3)+(k(1)*(i-k(2)))/(1+(i-k(2))^2);
    
    F(2) = F(2) + i*(k(3)+(k(1)*(i-k(2))))/(1+(i-k(2))^2);
    
    F(3) = F(3) + i^2*(k(3)+(k(1)*(i-k(2))))/(1+(i-k(2))^2);
    
end

F(2) = F(2) - mu_d;

F(3) = F(3) - var_d;

end