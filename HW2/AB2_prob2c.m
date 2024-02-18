%% Set simulation parameters
dt = 0.01; %Time-step
T = 10; %Final time
N = T/dt; %Number of steps to final time
%% Set ODE parameters
f = @(t,y) -y.^2; %The function on the rhs of the ODE
y0 = 1; %Initial value
y1 = 1/(1+dt); %Exact value at t_1

%[tvec,yvec] = AB2(0,y0,y1,f,dt,N);

% find 10 dt values from 10^-4 to 10^-1
dt_values = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5];
size(dt_values)
%initial the error varables
l2_errors_FE = zeros(size(dt_values));
l2_errors_RK2 = zeros(size(dt_values));
l2_errors_ab2 = zeros(size(dt_values));

%find errors in each dt values
for i = 1:length(dt_values)
    dt = dt_values(i);
    N = round(T/dt);
    y1 = 1/(1+dt);

    [tvec,yvec] = AB2(0,y0,y1,f,dt,N);
    %y_exact = 1./(1 + tvec)
    [tvec_FE,yvec_FE] = AB2_FE(0,y0,f,dt,N);
    [tvec_RK2,yvec_RK2] = AB2_RK2(0,y0,f,dt,N);
    
    %make the y can be compared at same points on the exact solution ab2.
    yvec_FE_compare = interp1(linspace(0, T, N+1), yvec_FE, tvec);
    yvec_RK2_compare = interp1(linspace(0, T, N+1), yvec_RK2, tvec);
    
    %calculating l2error for each method
    l2_errors_FE(i) = sqrt(sum((yvec - yvec_FE_compare).^2) / sum(yvec.^2));
    l2_errors_RK2(i) = sqrt(sum((yvec - yvec_RK2_compare).^2) / sum(yvec.^2));
   
end

loglog(dt_values, l2_errors_FE, '-o', dt_values, l2_errors_RK2, '-*');
title('FE-AB2 vs. RK2-AB2');
xlabel('\Delta t');
ylabel('Relative l2 Error');
legend('FE-AB2', 'RK2-AB2');
grid on;