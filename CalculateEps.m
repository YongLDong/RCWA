function [Rho,epsilon]=CalculateEps(Lambda,IncAngle,Psi,Delta,Title)
%calculate the epsilon and refractive index from ellipsoetry data
%unit of wavelength: Lambda [m]
%unit of incident angle: IncAngle [degree]

Rho=tand(Psi).*exp(-1i.*Delta.*pi./180);

epsilon=sind(IncAngle).^2+sind(IncAngle).^2.*tand(IncAngle).^2.*((1-Rho)./(1+Rho)).^2;

%%plot results
%dielectric constant
figure(1)
plot(Lambda.*1e9,real(epsilon))
hold on
plot(Lambda.*1e9,imag(epsilon))
legend('real part','imaginary part')
legend('boxoff')
FormatPlot('wavlength[nm]','dielectric constant',Title)

%pause for refractive index plot
figure(2)
RefIndex=sqrt(epsilon);
plot(Lambda.*1e9,real(RefIndex))
hold on
plot(Lambda.*1e9,imag(RefIndex))
legend('real part','imaginary part')
legend('boxoff')
FormatPlot('wavlength[nm]','refractive index',Title)
close
end