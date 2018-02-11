% In this script one solves the MDR (missing data recovery) problem
% A data vector f represented by a sparse (few nonzeros xs) expansion
% in a frame (here the discrete cosine transform DCT) Phi. f experiences
% a loss of information, its last n-m components are lost in b. The
% system Phi(1:m,:)*xr = b that would need to be solved to recover the
% coefficients xr is underdetermined and thus has infinitely many
% solutions. Since it is known that xs is sparse, one tries to
% find a sparse xr vector by solving the optimization problem
%
%     min ||x||_0  subject to   Ax = b, ||.||_0 denotes no. of nonzeros
%
% This is a combinatorial problem and NP hard. One replaces the ||.||_0
% norm by the ||.||_1 norm which encourages sparsity when minimized.
%
%     (MDR)     min ||x||_1  subject to   Ax = b  (x free)
%
% This problem can be rewritten as standard form LP and solved by
% Mehrotra's interior point method. Sparsity should be exploited
% Use m = 400, 350, 300
%
n = 500;			% data size
m = 350;			% available data size; may be changed
Phi = dctmtx(n);		% get the frame
xs = zeros(n,1);		% initialize x*
k = 10;				% number of nonzeros in x*
rand('state',1011);		% for reproducibility
p = randperm(60);		% indices of nonzeros in x*
randn('state',11);		% for reproducibility
xs(p(1:k)) = randn(k,1);	% random values for x*
f = Phi*xs;			% complete data vector v*
b = f(1:m);			% available data b
A = Phi(1:m,:);			% matrix A
%%
subplot(4,1,1);			% put all graphs in one figure
plot(1:n,f);			% plot complete data
subplot(4,1,2);
plot(1:n,[b;zeros(n-m,1)]);	% plot available data
%
% fill in here to solve for xr in Matlab or import your solution
% found with your program in a different language
%
load('x_opt.mat');
subplot(4,1,3);
plot(1:n,Phi*xr);		% plot recovered data
subplot(4,1,4);
plot(1:n,Phi*xr-f);		% plot error
