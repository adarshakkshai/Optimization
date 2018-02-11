minimize 
 obj: x1 + 2 x5 - x8
subject to
 c1:   3 x1 +   x2         -  2 x4    - x5                 -    x8 >= 2.5
 c2:          2 x2 + 1.1 x3                                        <=  2.1
 c3:                     x3                +  x6                    =  4.0
 c4:                        2.8 x4               -1.2 x7           >=  1.8
 c5:                        2.8 x4               -1.2 x7           <=  5.0
 c6: 5.6 x1                           + x5                + 1.9 x8 >=  3.0
 c7: 5.6 x1                           + x5                + 1.9 x8 <= 15.0
bounds  \default bounds are 0 <= x < +inf
 2.5 <= x1 < +inf
   0 <= x2 <= 4.1
 0.5 <= x5 <= 4.0
   0 <= x6 < +inf
   0 <= x7 < +inf  
   0 <=	x8 <= 4.3
end


\This is where I am performing branch and bound
\level 2 branch and bound
\x2 <= 1    	\x2 = 1.05 so I am bounding above 1
\2 <= x2 		\INFEASIBLE

\level 3 branch and bound
\x4 <= 0		\x4 = .643 so I am bounding above by 0, INFEASIBLE
\1 <= x4		

\level 3 branch and bound
\x5 <= 0		\x5 = 0.5, bounding above by 0 (NOT ALLOWED PER ORIGINAL CONSTRAINTS)
\1 <= x5			

\level 4 branch and bound
\x7 <= 0		\INTEGER SOLUTION!!!
\1 <= x7	

\level 5 branch and bound
\x4 <= 1		\INFEASIBLE
\2 <= x4			\NOT ALLOWED PER !!! (gives a solution ---> x7 fractional)

\level 6 branch and bound
\x7 <= 3		\INTEGER SOLUTION!!!
\4 <= x7			\fractional solution, now x2 and x4 are fractions

\level 7 branch and bound
\x2 <= 0		\INFEASIBLE
\1 <= x2			\combined with level 2 this implies x2 = 1
\x4 <= 2		\INFEASIBLE 
\3 <= x4			\INFEASIBLE

\level 8
\x4 <= 2		\INFEASIBLE
\3 <= x4		\INFEASIBLE
general
 x2 x3 x4 x5 x6 x7 x8
\binary
\ x3 x4
end
