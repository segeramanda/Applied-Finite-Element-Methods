function [A,M,Bvec] = assemble(p,t)
N = size(p,2);
A = sparse(N,N);
M = sparse(N,N); % allocate mass matrix
Bvec = zeros(N,1);

% assemble stiffness matrix A, mass matrix M and load vector b.
for K = 1:size(t,2)
    nodes = t(1:3,K);
    x1 = p(1,nodes);
    x2 = p(2,nodes);
 
    area_K = polyarea(x1,x2);
    b=[x2(2)-x2(3); x2(3)-x2(1); x2(1)-x2(2)]/2/area_K;
    c=[x1(3)-x1(2); x1(1)-x1(3); x1(2)-x1(1)]/2/area_K;
    
    AK = (b*b'+c*c')*area_K;                % element stiffness matrix
    A(nodes,nodes) = A(nodes,nodes) + AK;
    
    MK = [2 1 1;
        1 2 1;
        1 1 2]/12*area_K;                   % element mass matrix
    M(nodes,nodes) = M(nodes,nodes) + MK; % add element masses to M
    
        bK = [func(x1(1),x2(1));
        func(x1(2),x2(2));
        func(x1(3),x2(3))]/3*area_K;
    
    Bvec(nodes) = Bvec(nodes) + bK;
end
end
