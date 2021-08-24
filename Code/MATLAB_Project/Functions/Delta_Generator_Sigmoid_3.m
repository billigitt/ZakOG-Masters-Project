function [K, Delta, F_final] = Delta_Generator_Sigmoid_3(N, mu_d, var_d)

fun = @(x) f3d(x, N, mu_d, var_d);

K_0 = [0.5 N/2 1];

options = optimset('TolFun',1e-10,'MaxFunEvals',1e5,'Maxiter',1e5);

K = fsolve(fun, K_0, options);

Delta = Delta_From_K(K, N);

F_final = f3d(K, N, mu_d, var_d);

function [F] = f3d(k, N, mu_d, var_d)

F = zeros(1, 3);

i = 1:N;
   
F(1) = sum(k(1)*tanh(i-k(2)).*exp(-(i-k(2)).^2).*(i-k(3)).^4);

F(2) = sum(i.*(k(1)*tanh(i-k(2)).*exp(-(i-k(2)).^2).*(i-k(3)).^4))-mu_d;

F(3) = sum(i.^2.*(k(1)*tanh(i-k(2)).*exp(-(i-k(2)).^2).*(i-k(3)).^4))-var_d;

end

    function Delta = Delta_From_K(y, n)
        
        Delta = zeros(1, n);
        
        for i = 1:n
            
            Delta(i) = (y(1)*tanh(i-y(2)).*exp(-(i-y(2)).^2).*(i-y(3)).^4);
            
        end
        
    end

end