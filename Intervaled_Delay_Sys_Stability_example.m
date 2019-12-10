clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('..\..\Toolkit'));

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = {'sedumi'; 'gurobi'; 'cplex'; 'scip'};

%% Using example 1 of the referenced paper
% (LMI approach to exponential stability of linear systems with interval 
% time-varying delays)

A = [-0.075, 0.004; 0.005, -0.006];
D = [-0.065, 0.004; 0.003, -0.005];

% Interval bounds on delay:
h1 = 0.1; h2 = 0.5;

n = size(A,1);

%% Solving the LMI

a_l = 0.1;
a_u = 0.9;

prob_ever_feasible = 1;

for itr = 1:100
    a = (a_u + a_l) / 2;
    
    P = sdpvar(n,n);
    Q = sdpvar(n,n);
    R = sdpvar(n,n);
    U = sdpvar(n,n);
    
    S1 = sdpvar(n,n,'full');
    S2 = sdpvar(n,n,'full');
    S3 = sdpvar(n,n,'full');
    S4 = sdpvar(n,n,'full');
    S5 = sdpvar(n,n,'full');
    
    M_11 = (A')*P + P*A + 2*a*P - (exp(-2*a*h1)+exp(-2*a*h2))*R ...
           +0.5*S1*(eye(n)-A) + 0.5*(eye(n)-A')*(S1') + 2*Q;
    M_12 = exp(-2*a*h1)*R - S2*A; M_13 = exp(-2*a*h2)*R - S3*A;
    M_14 = P*D - S1*D - S4*A; M_15 = S1 - S5*A;
    M_22 = -exp(-2*a*h1)*(Q+R); M_24 = S2*D + exp(-2*a*h2)*U;
    M_33 = -exp(-2*a*h1)*(Q+R+U); M_34 = -S3*D + exp(-2*a*h2)*U;
    M_44 = 0.5*(S4*D + (S4*D)') - exp(-2*a*h2)*U;
    M_55 = S5 + S5' + (h1^2+h2^2)*R + ((h2-h1)^2)*U;
    
    constr = [ M_11,     M_12,     M_13,      M_14,    M_15;
              M_12',     M_22, zeros(n),      M_24,      S2;
              M_13', zeros(n),     M_33,      M_34,      S3;
              M_14',    M_24',    M_34',      M_44, S4-S5*D;
              M_15',      S2',      S3', (S4-S5*D),    M_55 ] <= 0;
          
    posit_constr = [P >= 0.1*eye(n); 
                    Q >= 0.1*eye(n); 
                    R >= 0.1*eye(n); 
                    U >= 0.1*eye(n)];
                
    constr = posit_constr + constr;
          
    ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                                 'solver',use_solver{1});

    optim = optimize(constr, [], ops);
    if optim.problem ~= 0 && optim.problem ~= 4
        if optim.problem == 1
            % Infeasible problem:
            a_u = a;
        else
            error(['Error solving problem: ' optim.info ])
        end
    else
        % Feasible problem:
        a_l = a;
        prob_ever_feasible = 1;
    end
    
    if (a_u-a_l) <= 1e-3
        break;
    end
end

if prob_ever_feasible
        fprintf('exponentially stable with alpha:\n');
        a = value(a)
        P = value(P)
        Q = value(Q)
        R = value(R)
        U = value(U)
else
    fprintf('\tNot exponentially stable\n');
end