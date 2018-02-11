var x1, >= 0, <= 50;	#cows
var x2, >= 0, <= 300;	#sheep
var y11, >= 0; var y12, >= 0;
var y21, >= 0; var y22, >= 0;

var z11; var z12; var z21; var z22;

param a11 = 1;
param a12 = .2;
param a21 = 150;
param a22 = 25;

param ahat11 = .05 * 1;
param ahat12 = .05 * .2;
param ahat21 = .05 * 150;
param ahat22 = .05 * 25;

param omega1 = .01;
param omega2 = .01;


maximize
obj: 250*x1 + 45*x2;

subject to

c1: a11*x1 + a12*x2 + ahat11*y11 + ahat12*y12 + omega1 * sqrt((ahat11*z11)^2 + (ahat12*z12)^2) <= 72;
c2: a21*x1 + a22*x2 + ahat21*y21 + ahat22*y22 + omega2 * sqrt((ahat21*z21)^2 + (ahat22*z22)^2)   <= 10000;

#index of x matches the column index of z and y
c3: x1 - z11 >= -y11;
c4: x1 - z11 <= y11;
c5: x1 - z21 >= -y21;
c6: x1 - z21 <= y21;

c7: x2 - z12 >= -y12;
c8: x2 - z12 <= y12;
c9: x2 - z22 >= -y22;
c10: x2 - z22 <= y22;
