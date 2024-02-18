%% Set simulation parameters
dt = 0.01; %Time-step
T = 10; %Final time
N = T/dt; %Number of steps to final time
%% Set ODE parameters
f = @(t,y) -y.^2; %The function on the rhs of the ODE
y0 = 1; %Initial value
y1 = 1/(1+dt); %Exact value at t_1

%% Run a simulation
%[tvec,yvec] = AB2(0,y0,y1,f,dt,N);
%[tvec,yvec] = AB2_FE(0,y0,f,dt,N);
[tvec,yvec] = AB2_RK2(0,y0,f,dt,N);
plot(tvec,yvec)