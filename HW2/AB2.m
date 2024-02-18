function [tvec,yvec] = AB2(t0,y0,y1,f,h,N)
% [tvec,yvec] = AB2(t0,y0,y1,f,h,N)
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
for n=1:N-1
   fvalue1 = f(tvec(n),yvec(n));
   fvalue2 = f(tvec(n+1),yvec(n+1));
   yvec(n+2) = yvec(n+1)+h/2*(3*fvalue2-fvalue1);
end