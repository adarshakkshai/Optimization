clear all;
clc;


%% Latitude and longitude are given in the form DDD.MM where DDD are the degrees and MM the minutes.
% A positive latitude is assumed to be \North", negative latitude means \South". Positive longitude means \East", 
% negative latitude is assumed to be \West". 
% For example, the input oordinates for Augsburg are 48.23 and 10.53, meaning 48o23 North and 10o53 East.

x = [16.47       96.1
16.47       94.44
20.09       92.54
22.39       93.37
25.23       97.24
22.00       96.05
20.47       97.02
17.20       96.29
16.30       97.38
14.05       98.12
16.53       97.38
21.52       95.59
19.41       97.13
20.09       94.55];

x = x';
n = numel(x(1,:));
x_cor = x(1,:);
y_cor = x(2,:);

deg = floor(x_cor);
degmin = x_cor - deg;
lat = pi.*(deg+5.0*degmin/3.0)/180;

deg = floor(y_cor);
degmin = y_cor - deg;
long = pi*(deg+5.0*degmin/3.0)/180;

RRR = 6378.388;
d = zeros(n); 

%% Generate Windy matrix.

for i = 1:n
    for j = i+1:n
        q1 = cos(long(i) - long(j));
        q2 = cos(lat(i) - lat(j));
        q3 = cos(lat(i)+lat(j));
        if (long(i)>long(j)) % long(i)>long(j) means that the longitude of city 'i' is much greater than the longitude of city 'j'. So, if I move to along path from 'i' to 'j', the distance gets longer 10% against west wind.   
           d(i,j) = floor(1.1.*(RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0)); % In case that I go from 'i' to 'j', the distance gets longer by 10%.
           d(j,i) = floor(0.9.*(RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0)); % In other cases, the distance gets shorter by 1-%.  
        elseif (long(j)>long(i))
           d(i,j) = floor(0.9.*(RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0));
           d(j,i) = floor(1.1.*(RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0));    
        else (long(i)==long(j))
           d(i,j) = floor((RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0)); 
           d(j,i) = floor((RRR.*acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0)); 
        end
    end
end


d_min = min(d(d~=0)); % min value of D matrix 
d_max = max(max(d));  % max value of D_matrix 

%% Generate D_prime matrix

if(~isnan(d_max)&&~isnan(d_min)) % (1) condition (page 4 of ATSP paper) = negative infinity < d_min <= d_max < positive infinity

e=0.1; % small positive number

D_prime=zeros(n,n);
D_prime_bar = zeros(2*n,2*n); % New symmetric matrix.


for i=1:n
    for j=1:n
        if i==j
            D_prime(logical(eye(size(d))))=0;                 % if i=j, then 0    1st condition 
        elseif (4*d_min-3*d_max>0)&&(i~=j)                    % if [4d_min-3d_max]>0 and i is not equal to j, then d(i,j)    
                D_prime(i,j) = d(i,j); 
        else                                                  % otherwise, d(i,j)+[3d_max-4d_min+e].
            D_prime(i,j) = d(i,j) +floor((3*d_max-4*d_min+e));        
        end
    end
end

end


%% Generate D_prime_bar matrix 

D_prime_max = max(max(D_prime));

D_prime_bar = [ floor(2.*n.*D_prime_max.*ones(n))        D_prime'
                D_prime         floor(2.*n.*D_prime_max.*ones(n)) ]; 
            
