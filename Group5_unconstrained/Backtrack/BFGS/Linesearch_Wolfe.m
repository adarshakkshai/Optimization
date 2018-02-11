function [ a_star ] = Linesearch_Wolfe(func,x0,p)

a = 0;  %initial alpha 
c2 = 0.9; % 0<c1<c2<1
c1 = 10e-4;

zoom1 = 0;
zoom2 = 0;  %counter variables to see how many times we call zoom and in which spot

amax = 1; % alpha max
r = .000001;
aa = amax*r; % current alpha

[phi_0, phi_grad_0] =  feval(func,x0,p,a);  %returns f(x0) and gradf(x0)
phi_a = phi_0;
phi_prime0 = phi_grad_0'*p;
i=1;
while (1)  
    [phi_aa, phi_grad_aa] = feval(func,x0,p,aa);  % Evaluate phi(alphai) and gradient of phi(alphai)
    phi_prime_aa = phi_grad_aa'*p;
    
    if (phi_aa>phi_0+c1*aa*phi_prime0) || (phi_aa>=phi_a && i>1)
        zoom1 = zoom1 + 1
        a_star=zoom(func,x0,p,a,aa,phi_0,phi_prime0,c1,c2);
        break
    end
    
    if abs(phi_prime_aa)<=-c2*phi_prime0
        a_star=aa;
        break
    end
    
    if phi_prime_aa>=0
        zoom2 = zoom2 + 1
        a_star=zoom(func,x0,p,aa,a,phi_0,phi_prime0,c1,c2);
        break 
    end
    
    a=aa;
    aa = (amax-a)*r + a; % increase to alpha(i+1)
    % update
    i=i+1;
	phi_a = phi_aa;
  
end



end

