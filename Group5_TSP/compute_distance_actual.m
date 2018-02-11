%this code takes the output vector from the fortran code main.f and
%newcdt.f that use INP.DAT to solve the assymetric TSP problem and then
%uses the assymetric distance matrix, d, to add up the "total distance"
%traveled in the tour
 
tour_vector = [8 14 4 5 6 12 13 2 10 1 9 7 11 3];
tour_distance = 0;

n = numel(tour_vector);
for i=1:n-1
    tour_distance = tour_distance + dist(min(tour_vector(i),tour_vector(i+1)), max(tour_vector(i),tour_vector(i+1)));
end
