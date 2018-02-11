n = 500;			% data size
m = 400;			% available data size; may be changed
Phi = dctmtx(n);		% get the frame
xs = zeros(n,1);		% initialize x*
k = 10;				% number of nonzeros in x*
rand('state',1011);		% for reproducibility
p = randperm(60);		% indices of nonzeros in x*
randn('state',11);		% for reproducibility
xs(p(1:k)) = randn(k,1);	% random values for x*
f = Phi*xs;			% complete data vector v*
b = f(1:m);			% available data b
A = Phi(1:m,:); % matrix A
%savefile = 'compression_sensing_data';
%save(savefile,'A','b','-ascii')
%A_writefile = 'A_compression_sensing_data_csv.txt'
%csvwrite(A_writefile,A);
%b_writefile = 'b_compression_sensing_data.txt';
%save(b_writefile,'b','-ascii')
%n_writefile = 'n.txt';
%save(n_writefile,'n','-ascii')
save('compression_data_n500_m400.mat','n','m','A','b');
save('compression_parameters_n500_m400.mat');
