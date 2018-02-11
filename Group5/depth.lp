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

\Adjust bounds to perform a depth first integer solution
\INFEASIBLE
\x2 = 1
\x3 = 0
\x4 = 1
\x5 = 1
\x6 = 4
\x7 = 0
\x8 = 1

\This yields an integer solution with obj = 4.5
x2 = 0
x3 = 0
x4 = 1 
x5 = 1
x6 = 4
x7 = 0
x8 = 0
end
