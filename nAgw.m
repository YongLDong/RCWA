function n = nAgw(wl)
%unit: wl[nm]
%%reading the real and imaginary part of refractive index (N=n+i*k)
%read data base from refractiveindex.info
[n_Cu,k_Cu,wls]=ReadRefractiveIndexCSV('./data/Ag.csv',2,1,50,2);
%interpolation of data
Wavelength=(0.21:0.0001:1.9);
n_Cu=interp1(wls,n_Cu,Wavelength);
k_Cu=interp1(wls,k_Cu,Wavelength);

Index=find(abs(Wavelength-wl/1e3)<1e-10);
n=n_Cu(Index)+1i.*k_Cu(Index);
end