function [f_grad] = extended_rosen_gradient(x)

alpha = 100;
n=length(x);
f_grad=zeros(n,1);

for i = 1:n/2 
    odd = -4*alpha*x(2*i-1)*(x(2*i)-x(2*i-1)^2)+2*(x(2*i-1)-1);
    even = 2*alpha*(x(2*i)-x(2*i-1)^2);
    f_grad(2*i-1) = odd;
    f_grad(2*i) = even;
end

end

 