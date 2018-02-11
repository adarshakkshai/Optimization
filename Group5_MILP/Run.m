clear all;
close all;
clc;

%% MIP Problem 
% 1. the number of iterations = 10 <--- Max objective val = -4.5
% 2. the number of iterations = 30 <--- Min objective val = 4.5

f = [-1, 0, 0, 0, -2, 0, 0, 1];                % f = objective function.
A = [-3, -1, 0, 2, 1, 0, 0, 1;                % A = inequality constraints 
     0, 2, 1.1, 0, 0, 0, 0, 0;
     0, 0, 0, -2.8, 0, 0, 1.2, 0;
     0, 0, 0, 2.8, 0, 0, -1.2, 0;
    -5.6, 0, 0, 0, -1, 0, 0, -1.9;
     5.6, 0, 0, 0, 1, 0, 0, 1.9];
b = [-2.5, 2.1, -1.8, 5, -3, 15];      % B = inequality constraints
Aeq = [0, 0, 1, 0, 0, 1, 0, 0];                % Aeq = equality constraints 
beq = 4;      % beq = equality constraints
lb = [2.5, 0, 0, 0, 0.5, 0, 0, 0];            % lb = lower bound
ub = [inf, 4.1, inf, inf, 4, inf, inf, 4.3];  % ub = upper bound 
M = [2 3 4 5 6 7 8]; 
e=1e-05;

[x_2 v_2 s_2] = IP1(f,A,b,Aeq,beq,lb,ub,M,e)