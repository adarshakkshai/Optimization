var x, >= 0, <= 50;
var y, >= 0, <= 300;

maximize
obj: 250*x + 45*y;

subject to

c1: x + .2*y <= 72;
c2: 150*x + 25*y <= 10000;
