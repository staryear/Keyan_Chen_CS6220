function [tvec,yvec] = AB3(t0,y0,y1,f,h,N)
% [tvec,yvec] = AB3(t0,y0,y1,f,h,N)
% Adams-Bashforth 2nd-order method
% Inputs
% t0,y0: initial condition
% y1:    additional start-up value
% f:     names of the right-hand side function f(t,y)
% h:     stepsize
% N:     number of steps
% Outputs
% tvec: vector of t values
% yvec: vector of corresponding y values
yvec = zeros(N+1,1);
tvec = linspace(t0,t0+N*h,N+1)';
yvec(1) = y0;
yvec(2) = y1;

% using RK3 to get third y value that will use in AB3
t1 = t0+h
k1 = f(t1, y1);
k2 = f(t1 + h/2, y1 + (1/2)*h*k1);
k3 = f(t1 + h, y1 - h*k1 + 2*h*k2);
yvec(3) = yvec(2) + (h/6)*(k1 + 4*k2 + k3);

for n=2:N-1
   %using first three value in y vector to calculate three fvalue in AB3
   fvalue1 = f(tvec(n-1), yvec(n-1));
   fvalue2 = f(tvec(n), yvec(n));
   fvalue3 = f(tvec(n+1), yvec(n+1));
   %using three f values to calculate next y value by using AB3 formula
   yvec(n+2) = yvec(n+1) + (h/12)*(23*fvalue3 - 16*fvalue2 + 5*fvalue1);
end