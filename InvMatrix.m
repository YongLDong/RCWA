function AInv=InvMatrix(A)
% Goal: solve A*x == b for x
% Set up some matrix A (I used a sparse matrix) -- do yourself
% Set up the vector b  -- do yourself
% Perform SVD on A
[U,S,V] = svd(A);
% A == U*S*V'  % Not needed, but you can check it yourself to confirm
% Calc number of singular values
s = diag(S);   % vector of singular values
% tolerance = max(size(A))*eps(max(s));
tolerance = 1e-14;
p = sum(s>tolerance);
% Define spaces
Up = U(:,1:p);
%U0 = U(:,p+1:Nx);
Vp = V(:,1:p);
%V0 = V(:,p+1:Nx);
%Sp = spdiags( s(1:p), 0, p, p );
SpInv = spdiags( 1.0./s(1:p), 0, p, p );
% Calc AInv such that x = AInv * b
AInv = Vp * SpInv * Up';
end