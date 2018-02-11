function [f, f_grad] = extended_rosen(x0,p,a)

n = length(x0);
x = x0 + a*p;
alpha = 100;
sum = 0;
for ii = 1:n/2
	xi = x(2*ii-1);
	xnext = x(2*ii);
	new = alpha*(xnext-xi^2)^2 + (xi-1)^2;
	sum = sum + new;
end

f = sum;

[f_grad] = extended_rosen_gradient(x);


end
