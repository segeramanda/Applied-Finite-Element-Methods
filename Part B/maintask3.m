clear all;
T = 2;
h = [1/5 1/20 1/40];
alpha = 4;
delta = 0.01;
%tvec = linspace(0,T,1000);
%k = tvec(2)-tvec(1); %time-step
geometry = @circleg;

for i = 1:length(h)
        
    hmax = h(i);
    timestep = 0.1*hmax;
    tvec = linspace(0,T, round(T/timestep));
    k = tvec(2)-tvec(1);
    
    
    [p,e,t] = initmesh(geometry ,'hmax',hmax);
    
    [A,M,b] = assemble(p,t);   %stiffness matrix A and mass matrix M
   
    for j = 1:length(p)
        u_0(j,1) = 1+20*rand();
    end
    u_former = u_0;
  
    for l = 1:length(tvec)-1
        LHS = (M+(delta/2)*k.*A);
        RHS = (M-(delta/2)*k.*A)*u_former(:,l)-k*M*(u_former(:,l)./(u_former(:,l)+alpha)-u_former(:,l)+(u_former(:,l)).^2);
        u_h = LHS\(RHS);             
        u_former(:,l+1) = u_h;
         
   % Population rates for every timestep  
        for K = 1:size(t,2)
            nodes = t(1:3,K);
            x = p(1,nodes);
            y = p(2,nodes);
            area_K = polyarea(x,y);
            Mp(K) = (u_h(nodes(1))+ u_h(nodes(2))... %Trapezoidal rule
                 +u_h(nodes(3)))/3*area_K; 
        end
        Mprey(l) = sum(Mp);
    end
figure()
%pdesurf(p,t,u_h)
pdeplot(p,[],t,'XYData',u_h,'XYStyle','interp',...
         'ZData',u_h,'ZStyle','continuous',...
         'ColorBar','on', 'ColorMap','default');
     
caption =sprintf('Solution for hmax=%.4f', hmax);
title(caption, 'Fontsize', 20)

figure(3)
plot(tvec(2:end),Mprey)
hold on;
caption =sprintf('Population rates');
legend('hmax=1/5', 'hmax = 1/20', 'hmax = 1/40')
title(caption, 'Fontsize', 20)
end

