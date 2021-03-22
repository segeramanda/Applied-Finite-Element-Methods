a = -1; % left end point of interval
b = 1; % right
N = 12; % number of intervals (start)
%h = (b-a)/N; % mesh size
%x = a:h:b; % node coords
x = linspace(a,b,N); % Vector of nodes
delta = 0.01;
TOL = 1e-3;
lambda = 0.9;
eta2 = 10; %start value to get in the while lop
while sum(sqrt(eta2))>TOL
    
    %Assemble matricies
    A = my_stiffness_matrix_assembler(x);
    B = my_load_vector_assembler(x,delta);
    M = my_mass_matrix_assembler(x);
    
    %Solving eq. systems
    xi = A\B;
    xi_biss = -M\(A*xi); %Discrete Laplcacian

    eta2 = my_eta2(x,xi_biss, delta); %error indicator 
    
    %Residuals
    R = zeros(length(x),1);
    for i = 1:length(x)
        R(i) = my_f(x(i)) + delta*xi_biss(i);
    end
    
    x1 =x;   %x value to plot with
    
    %Refinement
    for i = 1:length(eta2)
        if eta2(i) > lambda*max(eta2)
            x = [x (x(i+1)+x(i))/2];   
        end
    end
    x = sort(x);
    N=length(x);
end

%Plotting the results
subplot(2,2,1)
plot(x1,xi)
title ('Solution')
xlabel('x')
ylabel('u_h')

subplot(2,2,2)
plot(x1,sqrt(eta2))
title ('Error indicator')
xlabel('x')
ylabel('eta(u_h)')

subplot(2,2,3)
plot(x1(2:end),[1./diff(x1)])
title ('Mesh-size distribution')
xlabel('x')
ylabel('Distribution(x)')

subplot(2,2,4)
plot(x1,R)
title ('Residuals')
xlabel('x')
ylabel('R(u_h)')