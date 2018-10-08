function n = nAgw(wl)
%unit: wl[nm]
%%reading the real and imaginary part of refractive index (N=n+i*k)
%read data base from refractiveindex.info
[n_Ag,k_Ag,wls]=ReadRefractiveIndexCSV('./data/Ag.csv',2,1,50,2);
%interpolation of data
Wavelength=(0.21:0.0001:1.9);
n_Ag=interp1(wls,n_Ag,Wavelength);
k_Ag=interp1(wls,k_Ag,Wavelength);

Index=find(abs(Wavelength-wl/1e3)<1e-10);
n=n_Ag(Index)+1i.*k_Ag(Index);
end