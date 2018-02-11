function [a] = backtrack(func,x0,p)

%func is the function to be optmized
%x0 is the current point
%p is the search direction
%
%
%a is the step length that satisfies the algorithm

a = .01; %need a > 0
rho = .7;	%need 0<rho<1
c = .9;	% need 0 < c_lo < c < 1
[fk,fkgrad] = feval(func,x0,p,0);   %x0 = x0 +0*p
[fkk,fkkgrad] = feval(func,x0,p,a); %x = x0 + a*p

while(fkk > fk + c*a*fkgrad'*p)
	a = rho * a;
	fk = fkk;
    fkgrad = fkkgrad;
	[fkk,fkkgrad] = feval(func,x0,p,a);
end

