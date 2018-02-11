function f = testfun(x)

%Taken from problem 5) of professor Mittleman's global.pdf
%
%input x in [0,1]^n
%
%parameters
% zj in [0,1]^n, aj > 0, for j = 1,...,k
% 
% need aj >= 2 for f to be differentiable
%
% the structure of f implies that we can choose an values we want for f(zj), j = 1,...k  

%test case 1
%n = 3, k = 4
n=3;
k=4;
Z = [.1 .7 .6;
	.75 .14 .9;
	.6 .1 .7;
	1 0 .5];
fz = [4 -1 3 7]';

A = [2.5 1.5 4.1 2.7]';

fnum = 0;
p = 1;
for i = 1:k
	for j = 1:k
		if j ~= i
			p = p * norm(x-Z(j,:))^A(j);
		end
	end
	p = p * fz(i);
	fnum = fnum + p;
	p = 1;
end

fdenom = 0;
pp = 1;
for i = 1:k
	for j = 1:k
		if j ~= i
			pp = pp * norm(x-Z(j,:))^A(j);
		end
	end
	fdenom = fdenom + pp;
	pp = 1;
end

f = fnum / fdenom; %-fnum/fdenom for maximization

 
