% RCWA test 4
% for different azimuthal angle
% reproduce sheng, et. al. 1982
%
 clc; clear all;
tic();
% UNITS
meters = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
nanometers = 1e-9 * meters;
degrees = pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
lam0 = 647.1 * nanometers; %free space wavelength
theta = [0:1:70] * degrees;             % elevation angle
NN = length(theta);
phi   = 0 * degrees;             % azimuthal angle
pte   = 0;%1/sqrt(2);                       % amplitude of TE polarization
ptm   = 1;%1j/sqrt(2);                       % amplitude of TM polarization

% DEVICE PARAMETERS
ur1 = 1.0; %permeability in reflection region
er1 = 1; %permittivity in reflection region
ur2 = 1.0; %permeability in transmission region
er2 = -17.42+1i*0.58; %permittivity in transmission region

n_layers_unitcell = 2;
n_bilayers = 1;
urd = [1.0]; %permeability of device
erd = [-17.42+1i*0.58];

d = 60 * nanometers;
d2 = 50 * nanometers;

Lx = 1000 * nanometers; %period along x
Ly = 1000 * nanometers; %period along y

% RCWA PARAMETERS
Nx = 512; %number of point along x in real-space grid
Ny = round(Nx*Ly/Lx); %number of point along y in real-space grid
NH= 13;%number of harmonics
PQ = NH*[1 1]; %number of spatial harmonics along x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: BUILD DEVICE LAYERS ON HIGH RESOLUTION GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROSS SECTIONAL GRID
NL = 2;
% INITIALIZE LAYERS TO er AND ur
UR(1:Nx,1:Ny,1) = ur1;% zeros(Nx,Ny,NL);
ER(1:Nx,1:Ny,1) = er1;% zeros(Nx,Ny,NL);
UR(1:Nx,1:Ny,2) = urd(1);% zeros(Nx,Ny,NL);
ER(1:Nx,1:Ny,2) = erd(1);% zeros(Nx,Ny,NL);

w = 0.34; %length of one side of triangle
% g_height= 70 * nanometers;
% grating vector along x
f = 1;
nx = round(f*w*Nx);
nx1 = floor((Nx - nx)/2);
nx2 = nx1 + nx;

ER(nx1+1:nx2,1:Ny,1) = erd(1);

% ny = round(f*w*Ny);
% ny1 = floor((Ny - ny)/2);
% ny2 = ny1 + ny;
% 
% ER(1:Nx,ny1+1:ny2,1) = erd(2);

%grating vectors along y
% f = 1;
% ny = round(f*w*Ny);
% ny1 = floor((Ny - ny)/2);
% ny2 = ny1 + ny;
% 
% ER(1:Nx,ny1+1:ny2,1) = erd(1);
% % 
% nx = round(f*w*Nx);
% nx1 = floor((Nx - nx)/2);
% nx2 = nx1 + nx;
% 
% ER(nx1+1:nx2,1:Ny,1) = erd(2);

L = [d d2];

% %visualize the pattern
% figure(1)
% pcolor(1:Nx,1:Ny,ER(:,:,3));shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: COMPUTE CONVOLUTION MATRICES OF EACH LAYER OF DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
URC = zeros(PQ(1)^2,PQ(2)^2,NL);
ERC = zeros(PQ(1)^2,PQ(2)^2,NL);
for i = 1:NL
    URC(:,:,i) = convmat(UR(:,:,i),PQ(1),PQ(2));
    ERC(:,:,i) = convmat(ER(:,:,i),PQ(1),PQ(2));
end
% URC(:,:,3) = convmat(UR(:,:,3),PQ(1),PQ(2));
% ERC(:,:,3) = convmat(ER(:,:,3),PQ(1),PQ(2));

REF(1:NN) = 0;
TRN(1:NN) = 0;

for ii=1:NN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4: COMPUTE WAVE VECTOR EXPANSIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = eye(PQ(1)^2,PQ(2)^2);
    Z = zeros(PQ(1)^2,PQ(2)^2);

    n1 = sqrt(er1);
    n2 = sqrt(er2);
    k0 = 2*pi/lam0;
    kinc = n1*[sin(theta(ii))*cos(phi); sin(theta(ii))*sin(phi); cos(theta(ii))];


    p = [-floor(PQ(1)/2):+floor(PQ(1)/2)];         % indices along x
    q = [-floor(PQ(2)/2):+floor(PQ(2)/2)];         % indices along y

    Kx = kinc(1) - 2*pi*p/(k0*Lx);
    Ky = kinc(2) - 2*pi*q/(k0*Ly);

    [Ky, Kx] = meshgrid(Ky,Kx);

    Kx = diag(sparse(Kx(:)));
    Ky = diag(sparse(Ky(:)));

    Kzr = -(sqrt(ur1*conj(er1)*I-Kx^2-Ky^2));
    Kzt = conj(sqrt(ur2*conj(er2)*I-Kx^2-Ky^2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 5: COMPUTE EIGEN MODES OF FREE SPACE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Kz = conj(sqrte(I));
    Q = [Kx*Ky I+Ky^2; -(I+Kx^2) -Kx*Ky];
    W0 = [I Z; Z I];
    V0 = -1j*Q;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 6: INITIALIZE GLOBAL SCATTERING MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xx1 yy1] = size(V0);
    II = eye(xx1,yy1);
    ZZ = zeros(xx1,yy1);
    Sg.s11 = ZZ;
    Sg.s12 = II;
    Sg.s21 = II;
    Sg.s22 = ZZ;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 7: MAIN LOOP FOR S-MATRICES THROUGH LAYERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     S = scattering_matrix_ith(ERC(:,:,3),URC(:,:,3),Kx,Ky,k0,g_height,W0,V0);
%     Sg = Redheffer_star_product(Sg,S,II);
    for i = 1:NL
        S = scattering_matrix_ith(ERC(:,:,i),URC(:,:,i),Kx,Ky,k0,L(i),W0,V0);
        Sg = Redheffer_star_product(Sg,S,II);
    end
%     
%     for j = 1:n_bilayers-1
%         for i = 1:NL
%             S = scattering_matrix_ith(ERC(:,:,i),URC(:,:,i),Kx,Ky,k0,L(i),W0,V0);
%             Sg = Redheffer_star_product(Sg,S,II);
%         end
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 8: COMPUTE RELECTION SIDE CONNECTION S-MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Sr, Wref] = scattering_matrix_ref(er1,ur1,Kx,Ky,Kzr,Z,I,II,W0,V0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 9: COMPUTE TRANSMISSION SIDE CONNECTION S-MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [St, Wtrn] = scattering_matrix_trn(er2,ur2,Kx,Ky,Kzt,Z,I,II,W0,V0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 10: UPDATE GLOBAL SCATTERING MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sgg = Redheffer_star_product(Sr,Sg,II);
    Sgg = Redheffer_star_product(Sgg,St,II);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 11: COMPUTE REFLECTION AND TRANSMITTED FIELDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta = zeros(PQ(1)^2,1);
    delta(ceil(length(delta)/2),1) = 1;

    n_hat = [0; 0; 1];
    %%% cross product of k_inc and n_hat : k_inc x n_hat
    [k_inc_cross_n_hat, mag_k_inc_cross_n_hat] = vector_cross_product(kinc/k0,n_hat);

    %%% direction of TE polarized light
    if (theta(ii) == 0)
        a_hat_te = [0; 1; 0];
    else
        a_hat_te = [k_inc_cross_n_hat.x/mag_k_inc_cross_n_hat;
                    k_inc_cross_n_hat.y/mag_k_inc_cross_n_hat;
                    k_inc_cross_n_hat.z/mag_k_inc_cross_n_hat];
    end

    %%% cross product of a_hat_te and k_inc : a_hat_te x k_inc
    [a_hat_te_cross_k_inc, mag_a_hat_te_cross_k_inc] = vector_cross_product(a_hat_te,kinc/k0);

    %%% direction of TM polarized light
    a_hat_tm = [a_hat_te_cross_k_inc.x/mag_a_hat_te_cross_k_inc;
                a_hat_te_cross_k_inc.y/mag_a_hat_te_cross_k_inc;
                a_hat_te_cross_k_inc.z/mag_a_hat_te_cross_k_inc];

    %%% Polarization vector
    EP = pte*a_hat_te + ptm*a_hat_tm;

    esrc = zeros(2*PQ(1)^2,1);
    esrc(1:length(esrc)/2,1)=EP(1)*delta;
    esrc(length(esrc)/2+1:end,1)=EP(2)*delta;

    csrc = Wref\esrc;

    cref = Sgg.s11*csrc;
    ctrn = Sgg.s21*csrc;

    eref = Wref*cref;
    etrn = Wtrn*ctrn;

    rx = eref(1:length(eref)/2);
    ry = eref(length(eref)/2+1:end);
    rz = -Kzr\(Kx*rx+Ky*ry);

    tx = etrn(1:length(etrn)/2);
    ty = etrn(length(etrn)/2+1:end);
    tz = -Kzt\(Kx*tx+Ky*ty);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 12: COMPUTE DIFFRACTION EFFICIENCIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r2 = abs(rx).^2+abs(ry).^2+abs(rz).^2;
    R = real(-Kzr/kinc(3))*r2;
    
    r_s=real(-Kzr/kinc(3))*ry;
    r_s_matrix(ii)=sum(r_s);
    r_p=real(-Kzr/kinc(3))*sqrt(rx.^2+rz.^2);
    r_p_matrix(ii)=sum(r_p);
    
    if real(sum(rz))<0
        r_p_matrix(ii)=-(r_p_matrix(ii));
    end
    
    Rho(ii)=r_p_matrix(ii)./r_s_matrix(ii);
    
    REF_all(ii) = sum(R);
    REF_zero(ii) = R(floor(NH^2/2)+1);

    t2 = abs(tx).^2+abs(ty).^2+abs(tz).^2;
    T = real((ur1/ur2)*Kzt/kinc(3))*t2;
    TRN(ii) = sum(T);
end
toc()

Tan_Psi=abs(Rho);
Exp_Delta=Rho./abs(Rho);
Psi=atand(Tan_Psi);
Delta=atan2d(imag(Exp_Delta),real(Exp_Delta));
for i=1:length(Delta)
   if Delta(i)<-90
       Delta(i)=Delta(i)+360;
   end
end
Delta=180-Delta;
%plot ellipsometry results
plot(theta/degrees,Psi)
hold on
plot(theta/degrees,real(Delta))
FormatPlot('incident angle','','ellipsometry data')

% plot reflectance and transmittance
plot(theta/degrees,abs(r_p_matrix),theta/degrees,abs(r_s_matrix));
FormatPlot('incident angle','Reflection coefficient','')

close

[Rho,epsilon]=CalculateEpsFixedLam(lam0,theta/degrees,Psi,Delta,'CuSi grating surface');

