%xr = [-1.0749e-06; 0; 1.1515; 4.073e-6; 5.5565e-1;...
%    0; -4.617e-8; 0;2.519e-6; -1.4245e-6];

subplot(4,1,1);			% put all graphs in one figure
plot(1:n,f);			% plot complete data
subplot(4,1,2);
plot(1:n,[b;zeros(n-m,1)]);	% plot available data
%
% fill in here to solve for xr in Matlab or import your solution
% found with your program in a different language
%
load('x_opt_n500_m400.mat');
subplot(4,1,3);
plot(1:n,Phi*xr);		% plot recovered data
subplot(4,1,4);
plot(1:n,Phi*xr-f);		% plot error