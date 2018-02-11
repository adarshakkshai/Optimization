function [ a_star ] = zoom(func,x0,p,a_low,a_hi,phi_0,phi_grad_0,c1,c2)


while (1)
%     if abs(a_low - a_hi) < 1e-5 
%         a_star = aj;
%         break
%     end
	aj = (a_low + a_hi) / 2; 	%bisection
   % x_c = x0 + alpha_j*p;   
    [phi_aj, phi_grad_aj] = feval(func,x0,p,aj);  % phi(aj) and gradient of phi(aj)
   % x_lo = x0 + alpha_p*p;
    [phi_a_low, ~] = feval(func,x0,p,a_low);  % phi(a_low)
    if (phi_aj>phi_0+c1*aj*phi_grad_0)||(phi_aj>=phi_a_low)
        a_hi = aj;
    else
        if abs((phi_grad_aj)'*p)<=-c2*phi_grad_0
            a_star = aj;
            break
        end
        if  (phi_grad_aj'*p)*(a_hi-a_low)>=0
            a_hi = a_low;
        end
        a_low = aj; 
    end  
end

end
