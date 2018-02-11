minimize 
 obj: 0.5 x1 + 0.6667 x2 - 1.8 x3 + 5 x4 - 0.5 x5
subject to
 c1: x1 	+ x2 	+ x3 	+ x4 	<= 3
 c2: x1 	- 2 x2			>= 2
 c3: x1				+ x4	>= 1
 c4:  		 x1	+ x5	= 3.5 
 c5: 			x3		+ x5 <= 10
bounds  \default bounds are 0 <= x < +inf
 0 <= x1 <= +inf
 0 <= x2 <= +inf
 0 <= x3 <= +inf
 0 <= x4 <= +inf
general 
  x1 x2 x4
end
