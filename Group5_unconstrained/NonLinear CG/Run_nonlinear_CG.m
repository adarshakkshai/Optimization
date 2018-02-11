%% Run Nonlinear Conjugate Gradient Method (extended Rosen)
clear all
clc

N = 500;
x0 = zeros(N,1);

xres = nonlin_CG_PR('extended_Powell','extended_Powell_gradient',x0)
% 
% %% Run Nonlinear Conjugate Gradient Method (extended Powell)
% clear all
% clc
% 
% N = 500;
% x0 = zeros(N,1);
% 
% xres = nonlin_CG_FR('extended_Powell','extende_Powell_gradient',x0);
