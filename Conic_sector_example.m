clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('..\..\Toolkit'));

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = {'sedumi'; 'gurobi'; 'cplex'; 'scip'};

%% Define the system matrices

A = [-0.5, 0, 1; 0 -2 0; 0 0 -5];
B = [5 3; 6 5; 7 6];
C = [2 -2 0; 0 0 1];

n = size(A,1); 
m = size(B,2);
p = size(C,1);

D = zeros(p,m);

%% Find Cone using bisection:

% Define constraints:

c_u = 1e50;
c_l = 0;

prob_ever_feasible = 0;

for itr = 1:100
    c = (c_l+c_u)/2;
    
    r2 = sdpvar(1,1);
    P = sdpvar(n,n);

    positive_constr = [P>=0; r2 >= 0; r2 <= 1e5];

    feas_constr = [P*A+(A')*P+(C')*C P*B-c*(C')+(C')*D;
               (P*B-c*(C')+(C')*D)' (D')*D-c*(D+D')+(c^2 - r2)*eye()] <= 0;

    tot_constr = positive_constr + feas_constr;

    ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                                 'solver',use_solver{1});

    optim = optimize(tot_constr, [], ops);
%     clc;
%     fprintf('\tProblem: ');
    if optim.problem ~= 0 && optim.problem ~= 4
        if optim.problem == 1
            % Infeasible problem:
            c_l = c;
            
%             fprintf('Infeasible\n');
        else
            error(['Error solving problem: ' optim.info ])
        end
    else
        % Feasible problem:
        c_u = c;
%         fprintf('Feasible\n');
        prob_ever_feasible = 1;
    end
    
    if (c_u-c_l) <= 1e-3
        break;
    end
end

clc;
if prob_ever_feasible
        fprintf('Cone sector (c,r):\n');
        c
        r2 = value(sqrt(r2))
else
    fprintf('\tNo Cone Found\n');
end
    