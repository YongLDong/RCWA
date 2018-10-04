function [Rho,epsilon]=CalculateEpsFixedLam(Lambda,IncAngle,Psi,Delta,Title)
%calculate the epsilon and refractive index from ellipsoetry data
%unit of wavelength: Lambda [m]
%unit of incident angle: IncAngle [degree]

Rho=tand(Psi).*exp(1i.*Delta.*pi./180);

epsilon=sind(IncAngle).^2+sind(IncAngle).^2.*tand(IncAngle).^2.*((1-Rho)./(1+Rho)).^2;

%%plot results
%dielectric constant
plot(IncAngle,real(epsilon))
hold on
plot(IncAngle,imag(epsilon))
legend('real part','imaginary part')
legend('boxoff')
FormatPlot('incident angle','dielectric constant',Title)
end