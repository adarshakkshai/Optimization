param n := 15;
param k := 8;
set S := 1 .. n+1;
var x {S} binary;

minimize obj: x[n+1];

subject to c1: 2*sum{i in 1 .. n} x[i] + x[n+1] = 2*k+1;

