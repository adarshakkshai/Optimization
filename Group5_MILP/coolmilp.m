clear all;
close all;
clc;

%% MIP Problem 
% 1. the number of iterations = 10 <--- Max objective val = -4.5
% 2. the number of iterations = 30 <--- Min objective val = 4.5

f = [1/2, 2/3, -1.8, 5, -0.5];                % f = objective function.
A = [1, 1, 1, 1 0;                % A = inequality constraints 
     -1, 2, 0, 0 0;
     -1, 0, 0, -1 0;
     0 0 1 0 1];
b = [3, -2, -1 4];      % B = inequality constraints
Aeq = [1, 0, 0, 0, 1];                % Aeq = equality constraints 
beq = 3.5;      % beq = equality constraints
lb = [0, 0, 0, 0, 0];            % lb = lower bound
ub = [inf, inf, inf, inf, inf];  % ub = upper bound 
M = [1,2,4]; 
e=1e-05;

[x_2 v_2 s_2] = IP1(f,A,b,Aeq,beq,lb,ub,M,e)
