minimize
 obj: x7

subject to
 2x1 + 2x2 + 2x3 + 2x4 + 2x5 + 2x6 + x7 = 7 \n=6,k=3 

bounds
 0< x1 < 1
 0< x2 < 1
 0< x3 < 1
 0< x4 < 1
 0< x5 < 1
 0< x6 < 1
 0< x7 < 1

\begin branch and bound
\x3 < 0
1<x3
\x4<0
1<x4
\x6<0
1<x6
x1<0
\1<x1	\INFEASIBLE
x5<0
x2<0
end


