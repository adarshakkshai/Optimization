var x1 := 0, >= 0;
var x2 := 0, >= 0;
var y1 := 0, >= 0;
var y2 := 0, >= 0;
#create the variables to create equalities in inner problem constraints
var b1 >= 0;
var b2 >= 0;

#variables for the lambdas in the lagrangian function
var lambda1>=0;
var lambda2>=0;


minimize obj: -x1^2 - 3*x2 - 4*y1 + y2^2;

subject to
c1: x1^2 + 2*x2 <= 4;

#These are the KKT conditions of the inner problem

#from the two components of the gradient of the lagrangian function
c2: 2*y1 - (-lambda1 * 2	+ lambda2 * 3  	) = 0;
c3: -5   - ( lambda1   		- lambda2 * 4	) = 0;

#from the inequality constraints of the inner problem
c4: x1^2 - 2*x1 + x2^2 - 2*y1 + y2 + 3 - b1 = 0;
c5: x2 + 3*y1 - 4*y2 - 4 - b2 = 0;

c6: lambda1*(x1^2 - 2*x1 + x2^2 - 2*y1 + y2 + 3) = 0;
c7: lambda2*(x2 + 3*y1 - 4*y2 - 4) = 0;

