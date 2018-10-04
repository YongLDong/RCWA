function n = nSiw(wl)
%unit: wl[nm]
%%reading the real and imaginary part of refractive index (N=n+i*k)
%read data base from refractiveindex.info
[n_Si,k_Si,wls]=ReadRefractiveIndexCSV('./data/Si.csv',2,1,47,2);
%interpolation of data
Wavelength=(0.21:0.0001:1);
n_Si=interp1(wls,n_Si,Wavelength);
k_Si=interp1(wls,k_Si,Wavelength);

Index=find(abs(Wavelength-wl/1e3)<1e-10);
n=n_Si(Index)+1i.*k_Si(Index);


end