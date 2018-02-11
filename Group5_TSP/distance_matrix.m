%To compute the distance matrix for symmetric TSP problem
num_cities = 14;
x = [16.47 16.47 20.09 22.39 25.23 22.00 20.47 17.20 16.30 14.05 16.53 ...
21.52 19.41 20.09]; %latitudes
y = [96.10 94.44 92.54 93.37 97.24 96.05 97.02 96.29 97.38 98.12 97.38 ...
95.59 97.13 94.55]; %longitudes
RRR = 6378.388;

degx = floor(x);
degy = floor(y);
minx = x - degx;
miny = y - degy;
latitude = pi * (degx + 5.0 * minx / 3.0) / 180.0;
longitude = pi * (degy + 5.0 * miny / 3.0) / 180.0;

dist = zeros(numel(x));

for i = 1:numel(x)
    for j = i+1:numel(x)
        q1 = cos( longitude(i) - longitude(j));
        q2 = cos( latitude(i) - latitude(j));
        q3 = cos( latitude(i) + latitude(j));   
        dist(i,j) = floor(RRR * acos( 0.5 *((1.0 + q1)*q2 - (1.0 - q1)*q3) ) + 1.0);
    end
end



        
    
