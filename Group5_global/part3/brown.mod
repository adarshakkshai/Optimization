
param n := 10000;
set S = 1 .. n;
set SS = 1 .. n-1;
var X {i in S}:= (-1)^i;


minimize obj: sum {i in SS} ((X[i]*X[i])^(X[i+1]^2 + 1) + (X[i+1]^2) ^ (X[i]^2 + 1));

   
