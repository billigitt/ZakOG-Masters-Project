% This function is relevant to Q3 (a). This function is identical to
% ODE_SEIR_model_q2, except that we now have another input parameter:
% mintime. This enables us to run simulations from time not equal to zero.
% This can then be used to contruct piece-wise models with different
% parameters during different time periods.

function [Classes] = ODE_SEIR_model_q3(para,ICs,mintime,maxtime)

%Set tolerance for ode45
opts = odeset('RelTol',1e-5);

%Run ODE using ODE45
[t, pop] = ode45(@diff_SEIR_model_q3, mintime:1:maxtime, [ICs.S ICs.E ICs.I ICs.H ICs.R], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E',pop(:,2), 'I',pop(:,3),'H', pop(:, 4), 'R',pop(:,5),'t',t);


%Diff equations
function dPop = diff_SEIR_model_q3(t,pop,para)

S=pop(1);
E=pop(2);
I=pop(3);
H=pop(4);
R=pop(5);

dS = -para.beta*S.*I/para.N;
dE = para.beta*S.*I/para.N-para.sigma*E;
dI = para.sigma*E-para.gamma*I;
dH = para.p_H*para.gamma*I;
dR = (1-para.p_H)*para.gamma*I;


dPop = [dS; dE; dI; dH; dR];

end

end