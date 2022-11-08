%% Adam Ignaciuk 
% Gaussian Beams 
clear
close all

%constants 
lambda = 0.6*1e-9;
lambda2 = 0.5;

z1=10; % distance from the waist

%constants required for Gaussian Beam 
A0 = 1; %Amplitude 
%z0 = 3000; %Distance 
Nx = 2048; %size of matrix 
Ny = 2048; 
dx = lambda2/2; %
dy = lambda2/2; %
d = 300; %size of aperture (side or 2*radius)
z0=d/2; % Rayleigh range
prop_dist = 50;
n0 = 1;

%Circle aperture+signal 

x=(-Nx/2+1:Nx/2)*dx; %x points 
y=(-Ny/2+1:Ny/2)*dy; %y points 
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
aperture = zeros(Ny,Nx);
aperture(R<d/2)=1;
uin1 = aperture;  %signal entering the setup 

gauss = GaussianBeam2D(A0,z1,z0,Nx,dx,lambda2);
uout_theory = GaussianBeam2D(A0,z1+prop_dist,z0,Nx,dx,lambda2);
uout_num = AS_propagate(gauss,z1+prop_dist,lambda2,n0,dx);
figure('Color','w');
subplot(1,3,1); plot(x,abs(uout_num(:,Nx/2+1)).^2);title("Gaussian beam numerical");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(1,3,2); plot(x,abs(uout_theory(:,Nx/2+1)).^2);title("Gaussian beam theory");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(1,3,3); plot(x,abs(uout_num(:,Nx/2+1)).^2-abs(uout_theory(:,Nx/2+1)).^2);title("Error Bar");xlabel('x [um]');ylabel('Intensity[a.u]');

figure('Color','w');
subplot(1,3,1); imagesc(x,y,abs(uout_num));title("Gaussian beam numerical");xlabel('x [um]');ylabel('y [um]'); colorbar;
subplot(1,3,2); imagesc(x,y,abs(uout_theory));title("Gaussian beam theory");xlabel('x [um]');ylabel('y [um]'); colorbar;
subplot(1,3,3); imagesc(x,y,abs(uout_theory)-abs(uout_num));title("Difference");xlabel('x [um]');ylabel('y [um]'); colorbar;
