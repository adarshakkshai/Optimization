### VARIABLES & PARAMETERS
param k:=3;   
param n:=12;

var x{1..n,1..k};
var alpha;

### OBJECTIVE FUNCTION
maximize obj: alpha;

subject to constr1{i in 1..n}: sum{m in 1..k}(x[i,m]*x[i,m]) = 1;
subject to constr2{i in 1..n,j in 1..n}: if (i<j) then sum{m in 1..k}((x[i,m]-x[j,m])*(x[i,m]-x[j,m]))>=alpha^2;

	
option knitro_options "feastol=0 opttol=0 ms_enable=1 ms_maxsolves=10";
solve;
display x,alpha;
printf ("20.12%f \n",alpha);
