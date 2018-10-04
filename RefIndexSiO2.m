function [n]=RefIndexSiO2(wls)
%Sellmeier equation of SiO2
%wavelength unit is um
B1=0.6961663;
B2=0.4079426;
B3=0.8974794;
C1=0.0684043^2;
C2=0.1162414^2;
C3=9.896161^2;
n=sqrt(1+(B1.*wls.^2)./(wls.^2-C1)+(B2.*wls.^2)./(wls.^2-C2)+(B3.*wls.^2)./(wls.^2-C3));
end