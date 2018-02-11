
set S ordered;
param n := card {S};

set SS := 1 .. (2**n - 2);

set POW {k in SS} := {i in S: (k div 2**(ord(i)-1)) mod 2 = 1};


set LINKS := {i in S, j in S: ord(i) < ord(j)};

param cost {LINKS} >= 0;
var X {LINKS} binary;

minimize TotCost: sum {(i,j) in LINKS} cost[i,j] * X[i,j];

subj to SubtourElim {k in SS : 3 <= card{POW[k]} <=n-1}:
 sum {i in POW[k], j in POW[k]: (i,j) in LINKS} X[i,j] <= card{POW[k]}-1;

subj to Tour: sum {(i,j) in LINKS} X[i,j] = n;


