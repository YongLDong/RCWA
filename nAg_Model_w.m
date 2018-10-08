function n = nAg_Model_w(wl)
%unit: wl[nm]
%%reading the real and imaginary part of refractive index (N=n+i*k)
%read data base from refractiveindex.info
Omg=2*pi./(wl*1e-3);
Eps_inf=1;
Delta_Eps_k=[1759.471 135.344 258.1946 22.90436 1749.06 11756.18];
A_k=[1 1 1 1 1 1];
B_k=[0.243097 19.68071 2.289161 0.329194 4.639097 12.25];
C_k=[0 17.07876 515.022 1718.357 2116.092 10559.42];
%interpolation of data
Eps=Eps_inf+zeros(1,length(wl));
for NN=1:6
    Eps=Eps+Delta_Eps_k(NN)./(A_k(NN)*(1i*Omg).^2-B_k(NN)*(1i*Omg)+C_k(NN));
end
n=sqrt(Eps);
end
% n=n_Ag(Index)+1i.*k_Ag(Index);