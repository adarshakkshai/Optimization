function [f_grad] = extended_Powell_gradient(x)

n=length(x);
f_grad=zeros(n,1);

for i = 1:n/4 % i=3
    x_1 = x(4*i-3);
	x_2 = x(4*i-2);
    x_3 = x(4*i-1);
    x_4 = x(4*i);
    
    Fir = 2*(x_1-10*x_2)+40*(x_1-x_4)^3;
    Sec = -20*(x_1-10*x_2)+4*(x_2-2*x_3)^3;
    Thi = 10*(x_3-x_4)-8*(x_2-2*x_3)^3;
    Fou = -10*(x_3-x_4)-40*(x_1-x_4)^3;
    
    f_grad(i) = Fir;
    f_grad(i+1) = Sec;
    f_grad(i+2) = Thi;
    f_grad(i+3) = Fou;
end

end

 