clc

clear all
close all

R = 2;

a = 1;
b = 5;

tau = 7;

total_time = 100;

change = 80;

I_0 = 3;

I = I_0;

w_o = [.1 .2 .3 .4];

w_a = [.5 .2 .2 .1];

w = [repmat(w_o, change, 1); repmat(w_a, total_time-change, 1)];

for i = 1:change-1

    I_new = poissrnd(R*Incidence_Generator_2(I, w_o));
    
    I = [I I_new];
    
end

for i = change:total_time
    
    I_new = poissrnd(R*Incidence_Generator_2(I, w_a));
    
    I = [I I_new];
    
end

Shape = zeros(total_time, 1);
Scale = zeros(total_time, 1);
Mean = zeros(total_time, 1);

for i = tau:total_time
     
    Shape(i) = a + sum(I(i-tau+1:i));
    
    for k = i:-1:i-tau+1
       
        I_relevant = I(1:k);
        
        Scale(i) = Scale(i) + Incidence_Generator_2(I_relevant, [0 w(k, :)]);
        
    end
    
    Scale(i) = 1/(Scale(i)+1/b);
    
    Mean(i) = Shape(i)*Scale(i);
    
end

figure(1)
plot(Mean)
figure(2)
plot(I)