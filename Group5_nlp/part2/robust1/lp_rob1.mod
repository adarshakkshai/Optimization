var x1, >= 0, <= 50;	#cows
var x2, >= 0, <= 300;	#sheep
var y1, >= 0;
var y2, >= 0;

param a11 = 1;
param a12 = .2;
param a21 = 150;
param a22 = 25;

param ahat11 = .25 * 1;
param ahat12 = .25 * .2;
param ahat21 = .25 * 150;
param ahat22 = .25 * 25;


maximize
obj: 250*x1 + 45*x2;

subject to

c1: a11*x1 + a12*x2 + ahat11*y1 + ahat12*y2 <= 72;
c2: a21*x1 + a22*x2 + ahat21*y1 + ahat22*y2 <= 10000;

c3: x1 >= -y1;
c4: x1 <= y1;
c5: x2 >= -y2;
c6: x2 <= y2;
