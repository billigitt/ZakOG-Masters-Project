function [K, Delta, F_final] = Delta_Generator_Sinusoidal_1(N, mu_d, var_d)

fun = @(x) f3d(x, N, mu_d, var_d);

K_0 = [1 1 1];

options = optimset('TolFun',1e-10,'MaxFunEvals',1e5,'Maxiter',1e5);

K = fsolve(fun, K_0, options);

Delta = Delta_From_K(K, N);

F_final = f3d(K, N, mu_d, var_d);

function [F] = f3d(k, N, mu_d, var_d)

F = zeros(1, 3);

i = 1:N;
   
F(1) = sum(k(1)*sin(2*pi*((i-k(3))/(N-k(3))).^2).*i.^2.*(i-k(2)).^2);

F(2) = sum(i.*(k(1)*sin(2*pi*((i-k(3))/(N-k(3))).^2).*i.^2.*(i-k(2)).^2))-mu_d;

F(3) = sum(i.^2.*(k(1)*sin(2*pi*((i-k(3))/(N-k(3))).^2).*i.^2.*(i-k(2)).^2))-var_d;

end

    function Delta = Delta_From_K(y, n)
        
        Delta = zeros(1, n);
        
        for i = 1:n
            
            Delta(i) = y(1)*sin(2*pi*((i-y(3))/(n-y(3))).^2).*i.^2.*(i-y(2)).^2;
            
        end
        
    end

end

