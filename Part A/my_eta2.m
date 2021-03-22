% error indicator "eta"
function eta2 = my_eta2(x, xi_biss, delta)
N = length(x);
eta2 = zeros(N,1); 
for i = 1:N-1
h = x(i+1)-x(i);

eta2(i) = eta2(i)+ 1/2*h^3*((my_f(x(i))+delta*xi_biss(i)).^2 +(my_f(x(i+1))+delta*xi_biss(i+1)).^2); 

end
end