clear all;
geometry = @circleg;
h = [1/2 1/4 1/8 1/16 1/32];
EnE = zeros(length(h),1);

for i = 1:length(h)
    hmax = h(i);
    [p,e,t] = initmesh(geometry ,'hmax',hmax);
    [A,M,b] = assemble(p,t);   %stiffness matrix
    
    u_e = @sin;
    x1 = p(1,:);
    x2 = p(2,:);
    u_exact = u_e(2*pi*x1).*u_e(2*pi*x2);
%     for j = 1:length(p)
%         x1 = p(1,j);
%         x2 = p(2,j);
%         u_exact(j,1) = u_e(2*pi*x1).*u_e(2*pi*x2);
%     end
    
    for k = 1:length(e(1,:))
        x11 = p(1,e(1,k));
        x22 = p(2,e(1,k));
        RHS(k,1) = u_e(2*pi*x11).*u_e(2*pi*x22);
    end
    
    I = eye(length(p));         %construct the identity matrix 
    A(e(1,:),:) = I(e(1,:),:);  %replace the rows corresponding 
                                    %to the boundary nodes by
                                        % corresponding rows of I
    b(e(1,:)) = RHS(:,1);       %put the boundary value into the RHS

    u_h = A\b;                  %solving linear eq.

    err=u_exact'-u_h;            %Error
    EnE(i,1)=sqrt(err'*A*err);  %Energy norm
    
    figure();
    pdeplot(p,[],t,'XYData',u_h,'XYStyle','interp',...
         'ZData',u_h,'ZStyle','continuous',...
         'ColorBar','on', 'ColorMap','default');
    %pdesurf(p,t,u_h)
    caption =sprintf('Solution for hmax=%.4f',hmax);
    title(caption, 'Fontsize', 20) 
end
figure(6)
    P = polyfit(log(h)',log(EnE),1);
    P = P(1);      %convergence rate
    
%     P = zeros(length(EnE)-2,1);
%     for n = 2:(length(EnE)-1)
%         P(n-1) = log(EnE(n+1)/EnE(n))/log(EnE(n)/EnE(n-1));       
%     end

loglog(h,EnE, 'b', h,h.^P, 'r')
title('Energy norm vesus hmax and hmax^p versus hmax, P = 1.4089')
legend('Ene', 'hmax^p')






