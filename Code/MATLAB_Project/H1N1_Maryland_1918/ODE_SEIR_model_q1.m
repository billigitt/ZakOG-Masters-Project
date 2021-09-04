function [Classes] = ODE_SEIR_model_q1(para,ICs,maxtime)

%Set tolerance for ode45, producing a 'finer' integration approximation
opts = odeset('RelTol',1e-5);

%Run ODE using ODE45, which calls the diff_SEIR_model_q1 function below as
%a handle. This is effectivelyu the integrand. The other input arguments
%are self-explanatory.
[t, pop] = ode45(@diff_SEIR_model_q1, 0:1:maxtime, [ICs.S ICs.E ICs.I ICs.R], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E',pop(:,2), 'I',pop(:,3),'R',pop(:,4),'t',t);


%Diff equations
function dPop = diff_SEIR_model_q1(t,pop,para)

S=pop(1);
E=pop(2);
I=pop(3);
R=pop(4);

dS = -para.beta*S.*I/para.N;
dE = para.beta*S.*I/para.N-para.sigma*E;
dI = para.sigma*E-para.gamma*I;
dR = para.gamma*I;


dPop = [dS; dE; dI; dR];

end

end