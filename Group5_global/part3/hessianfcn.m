function H = hessianfcn(x,~)

% Evaluate the (sparse, symmetric) Hessian matrix for the function
% defined in brownfgh.m

%This code was copied from the evaluate hessian section of brownfgh.m
n = length(x);
v=zeros(n,1);
i=1:(n-1);
v(i)=2*(x(i+1).^2+1).*((x(i).^2).^(x(i+1).^2))+...
    4*(x(i+1).^2+1).*(x(i+1).^2).*(x(i).^2).*((x(i).^2).^((x(i+1).^2)-1))+...
    2*((x(i+1).^2).^(x(i).^2+1)).*(log(x(i+1).^2));
v(i)=v(i)+4*(x(i).^2).*((x(i+1).^2).^(x(i).^2+1)).*((log(x(i+1).^2)).^2);
v(i+1)=v(i+1)+...
    2*(x(i).^2).^(x(i+1).^2+1).*(log(x(i).^2))+...
    4*(x(i+1).^2).*((x(i).^2).^(x(i+1).^2+1)).*((log(x(i).^2)).^2)+...
    2*(x(i).^2+1).*((x(i+1).^2).^(x(i).^2));
v(i+1)=v(i+1)+4*(x(i).^2+1).*(x(i+1).^2).*(x(i).^2).*((x(i+1).^2).^(x(i).^2-1));
v0=v;
v=zeros(n-1,1);
v(i)=4*x(i+1).*x(i).*((x(i).^2).^(x(i+1).^2))+...
    4*x(i+1).*(x(i+1).^2+1).*x(i).*((x(i).^2).^(x(i+1).^2)).*log(x(i).^2);
v(i)=v(i)+ 4*x(i+1).*x(i).*((x(i+1).^2).^(x(i).^2)).*log(x(i+1).^2);
v(i)=v(i)+4*x(i).*((x(i+1).^2).^(x(i).^2)).*x(i+1);
v1=v;
i=[(1:n)';(1:(n-1))'];
j=[(1:n)';(2:n)'];
s=[v0;2*v1];
H=sparse(i,j,s,n,n);
H=(H+H')/2;
end