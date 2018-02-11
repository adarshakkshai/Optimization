
set S ordered;
param n := card {S};

set LINKS := {i in S, j in S: ord(i) < ord(j)};

param cost {LINKS} >= 0;
var X {LINKS} binary;

minimize TotCost: sum {(i,j) in LINKS} cost[i,j] * X[i,j];

#this ensures that each node only has two distinct edges connected to it
subj to Tour {i in S}: 
   sum {(i,j) in LINKS} X[i,j] + sum {(j,i) in LINKS} X[j,i] = 2;

#Make sure that every node is used
subj to Tour1: sum {(i,j) in LINKS} X[i,j] = n;

#First set of subtours that we need to eliminate
subj to elim1: X['a','h'] + X['a','b'] + X['b','h'] <=2;

subj to elim2: X['c','d'] + X['d','e'] + X['e','f'] +
	X['f','l'] + X['g','l'] + X['g','m'] + X['m','n'] +
   	X['c','n'] <= 7;

subj to elim3: X['i','j'] + X['j','k'] + X['i','k'] <= 2;

#Second set of subtours that we need to eliminate
subj to elim4: X['a','b'] + X['b','j'] + X['i','j'] +
	X['i','k'] + X['h','k'] + X['a','h'] <= 5;

subj to elim5: X['c','d'] + X['d','e'] + X['e','l'] +
	X['f','l'] + X['f','g'] + X['g','m'] + X['m','n'] +
	X['c','n'] <= 7;

#Third set of subtours that we need to eliminate
subj to elim6: X['a','b'] + X['b','h'] + X['h','k'] +
	X['i','k'] + X['i','j'] + X['a','j'] <= 5;

subj to elim7: X['c','d'] + X['d','e'] + X['e','f'] +
	X['f','l'] + X['l','m'] + X['g','m'] + X['g','n'] +
	X['c','n'] <= 7;

#Fourth set of subtours that we need to eliminate
subj to elim8: X['a','b'] + X['b','j'] + X['j','k'] +
	X['i','k'] + X['h','i'] + X['a','h'] <= 5;

subj to elim9: X['c','d'] + X['d','l'] + X['f','l'] +
	X['e','f'] + X['e','g'] + X['g','m'] + X['m','n'] +
	X['c','n'] <= 7;

#Fifth set of subtours that we need to eliminate
subj to elim10: X['a','h'] + X['b','h'] + X['b','j'] +
	X['i','j'] + X['i','k'] + X['a','k'] <= 5;

subj to elim11: X['c','d'] + X['d','e'] + X['e','l'] +
	X['f','l'] + X['f','m'] + X['g','m'] + X['g','n'] +
	X['c','n'] <= 7;

#Sixth set of subtours that we need to eliminate
subj to elim12: X['a','b'] + X['b','h'] + X['h','j'] +
	X['i','j'] + X['i','k'] + X['a','k'] <= 5;

subj to elim13: X['c','d'] + X['d','f'] + X['f','l'] +
	X['e','l'] + X['e','g'] + X['g','m'] + X['m','n'] +
	X['c','n'] <= 7;





