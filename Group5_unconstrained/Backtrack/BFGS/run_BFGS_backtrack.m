N = 500;
x0 = zeros(N,1);
xres = BFGS_backtrack('extended_rosen',x0)   %if this works correctly then xres != ones(N,1)

