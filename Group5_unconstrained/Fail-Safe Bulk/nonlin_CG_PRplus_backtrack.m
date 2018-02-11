function [xres] = nonlin_CG_PRplus(func,funcgrad,x0)
% func evaluates both the value of the function and the gradient;
% func evaluates only the gradient of the function;

f0_grad = feval(funcgrad,x0); 
p0 = -f0_grad;
k = 0;
l2 = norm(f0_grad);
tol = 1e-5;


while abs(l2) > tol
        a = backtrack(func,x0,p0);   % we need to perform a line search in place of the formula (5.24a)
        xres = x0 + a*p0;
        f_grad=feval(funcgrad,xres);
        
        betaPR = (f_grad'*(f_grad-f0_grad)) / (f0_grad'*f0_grad);
        betaPR = max(0,betaPR);
		p1 = -f_grad+betaPR*p0;
         
        p0 = p1;
        k=k+1;
        f0_grad=f_grad;
        l2=norm(f_grad);
        x0=xres;
        k
end


%% while ---> for
% for k=0:400
%         a = Linesearch_Wolfe(functname,p0,x0);   % we need to perform a line search in place of the formula (5.24a)
%         xres = x0 + a*p0;
%         f_grad=extended_rosen_gradient(xres);
%         
%         beta = (f_grad'*f_grad)/(f0_grad'*f0_grad);
%         p1 = -f_grad+beta*p0;
%          
%         p0 = p1;
%         f0_grad=f_grad;
%         l2=norm(f_grad);
%         x0=xres;
% end



end

