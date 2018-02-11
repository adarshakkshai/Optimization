param n := 20;
var X {i in 1..n};
let X[1] := -1;
let X[2] := 1;
let X[3] := 1;
let X[4] := 1;
let X[5] := 1;
let X[6] := 1;
let X[7] := 1;
let X[8] := 1;
let X[9] := 1;
let X[10] := 1;
let X[11] := 1;
let X[12] := 1;
let X[13] := 1;
let X[14] := 1;
let X[15] := 1;
let X[16] := 1;
let X[17] := 1;
let X[18] := 1;
let X[19] := 1;
let X[20] := 1;

minimize obj: sum{i in 1..n-1} ( 100 * (X[i+1] - X[i]^2)^2 + (1-X[i])^2 );

