% This function is relevant to Q2 (a). Everything is identical to
% ODE_SEIR_model_q1.m, except that we introduce an H class into the
% diff_SEIR_model_q2 function. Technically this is an SEIRH model.

function [Classes] = ODE_SEIR_model_q2(para,ICs,maxtime)

%Set tolerance for ode45
opts = odeset('RelTol',1e-5);

%Run ODE using ODE45
[t, pop] = ode45(@diff_SEIR_model_q2, 0:1:maxtime, [ICs.S ICs.E ICs.I ICs.H ICs.R], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E',pop(:,2), 'I',pop(:,3),'H', pop(:, 4), 'R',pop(:,5),'t',t);


%Diff equations
function dPop = diff_SEIR_model_q2(t,pop,para)

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