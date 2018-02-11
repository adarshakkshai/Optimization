
%fmincon finds the contrained minimum of a function of several variables
%
%[x] = fmincon(@func, x0,A,B,Aeq,Beq,LB,UB,,NONLCON,options)
%
% @func is the function to be optimized
% x0 is the initial point
% A,B are inequality contraints AX <= B 
% Aeq, Beq are equality constraints AeqX = Beq
% LB and UB are lower and upper bounds

%specify interior-point algorithm (this is actually the default algorithm)
%We are not going to use the gradient or the hessian

options = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',[1000000]);

n = 10000;
x0 = ones(n,1);
x0(2:2:end) = -x0(2:2:end);
tic
[x,obj] = fmincon(@brownfgh, x0, [],[],[],[],[],[],[],options);
toc
obj
