function [f, f_grad] = extended_Powell(x,p,a)
% x=x0;

n = length(x);
x = x + a*p; %we need to evaluate with the given step direction and step length
sum = 0;
for i = 1:n/4
	x_1 = x(4*i-3);
	x_2 = x(4*i-2);
    x_3 = x(4*i-1);
    x_4 = x(4*i);
	new = (x_1-10*x_2)^2+5*(x_3-x_4)^2+(x_2-2*x_3)^4+10*(x_1-x_4)^4;
	sum = sum + new;
end

f = sum;

[f_grad] = extended_Powell_gradient(x);


end
