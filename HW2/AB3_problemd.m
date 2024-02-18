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
dt_values = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001];
%dt_values = [0.01]
size(dt_values)
%initial the error varables
l2_errors_AB3 = zeros(size(dt_values));
l2_errors_AB2 = zeros(size(dt_values));

%find errors in each dt values
for i = 1:length(dt_values)
    dt = dt_values(i);
    N = round(T/dt);
    y1 = 1/(1+dt);

    [tvec_AB2,yvec_AB2] = AB2(0,y0,y1,f,dt,N);
    [tvec_AB3,yvec_AB3] = AB3(0,y0,y1,f,dt,N);
    y_exact = 1./(1 + tvec_AB2)
    
    %make the y can be compared at same points on the exact solution ab2.
    %yvec_AB3_compare = interp1(linspace(0, T, N+1), yvec_AB3, tvec_AB2);
    %yvec_AB2_compare = interp1(linspace(0, T, N+1), yvec_AB2, tvec_AB2);
    
    %calculating l2error for each method
    l2_errors_AB3(i) = sqrt(sum((y_exact - yvec_AB3).^2) / sum(y_exact.^2));
    l2_errors_AB2(i) = sqrt(sum((y_exact - yvec_AB2).^2) / sum(y_exact.^2));
   
end

loglog(dt_values, l2_errors_AB2, '-o', dt_values, l2_errors_AB3, '-*');
%loglog(tvec, errors_AB2, '-o', tvec, errors_AB3, '-*');
title('AB2 vs. AB3');
xlabel('\Delta t');
ylabel('Relative l2 Error');
legend('AB2', 'AB3');
grid on;