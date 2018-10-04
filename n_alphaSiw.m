function n = n_alphaSiw(wl)
%unit: wl[nm]
%%reading the real and imaginary part of refractive index (N=n+i*k)
%read data base from refractiveindex.info
[n_alphaSi,k_alphaSi,wls]=ReadRefractiveIndexCSV('./data/alpha_Si.csv',2,1,47,2);
%interpolation of data
Wavelength=(0.21:0.0001:1);
n_alphaSi=interp1(wls,n_alphaSi,Wavelength);
k_alphaSi=interp1(wls,k_alphaSi,Wavelength);

Index=find(abs(Wavelength-wl/1e3)<1e-10);
n=n_alphaSi(Index)+1i.*k_alphaSi(Index);


end