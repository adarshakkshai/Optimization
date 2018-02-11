function xres = BFGS_backtrack(func,x0)

N = numel(x0);
I = eye(N,N);
Hk = I;
zero = zeros(N,1);
[~,z] = feval(func,x0,zero,0);	%return the function and the gradient evaluated at x0 
l2 = norm(z);
xprev = x0
tol = 1e-5;
k = 0;

%%
while  l2 > tol
	
    p = -Hk*z;
    a = backtrack(func,xprev,p);	
    xnext = xprev + a*p;
    s = xnext - xprev;
	[~,zz] = feval(func,xnext,zero,0);		%the gradient evaluated at x(k+1)
    y = zz - z;
    Hk = (I - (1/(y'*s))*s*y')*Hk*(I - (1/(y'*s)*y*s')) + (1/(y'*s))*(s*s');
    xprev = xnext;
	xres = xprev
    l2 = norm(zz);
    k = k+1;
	z = zz;
    k
end
