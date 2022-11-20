%% Error analysis 

clear
close all

%constants 
lambda = 0.6*1e-9;
lambda2 = 0.5;

z1=(0:0.5:15); % distance from the waist

%constants required for Gaussian Beam 
A0 = 1; %Amplitude 
%z0 = 3000; %Distance 
Nx = 2048; %size of matrix 
Ny = 2048; 
dx = lambda2/2; %
dy = lambda2/2; %
d = 1; %size of aperture (side or 2*radius)
z0=d/2; % Rayleigh range
prop_dist = 50;
n0 = 1;

x=(-Nx/2+1:Nx/2)*dx; %x points 
y=(-Ny/2+1:Ny/2)*dy; %y points 
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
aperture = zeros(Ny,Nx);
aperture(R<d/2)=1;
uin1 = aperture;  %signal entering the setup 


%uout1 = [];
%uout_theory = [];
uoutmiddle = [];
uout_theorymiddle= [];
gauss = GaussianBeam2D(A0,z1(1),z0,Nx,dx,lambda2);
for i = 0:0.5:15

    uout = AS_propagate(gauss, i+prop_dist, lambda2, n0, dx);
    uoutmiddle = [uoutmiddle uout(1024,1024)];
    
    uth = GaussianBeam2D(A0,i+prop_dist,z0,Nx,dx,lambda2);
    uout_theorymiddle = [uout_theorymiddle uth(1024,1024)];

end 

I_n = (abs(uoutmiddle)).^2;
I_t = (abs(uout_theorymiddle)).^2;

figure(Color='w');
subplot(1,2,1); plot(z1,I_t); xlabel("z[m]"); ylabel("I theoretical [a.u.]"); title("I of central point[theory]");
subplot(1,2,2); plot(z1,I_n); xlabel("z[m]"); ylabel("I numerical[a.u.]"); title("I of central point[numerical]")

sig = (abs(uout_theorymiddle - uoutmiddle)).^2;

figure(Color='w');
plot(z1,sig); xlabel("z[m]"); ylabel("Error"); title("Error of numerical propagation");

uanali1 = [];
uanali5 = [];
uanali75 = [];
z2 = (0:0.05:15)
for z = 0:0.05:15
    a1 = Analitical(lambda2,z,1);
    uanali1 = [uanali1 abs(a1).^2];
    a2 = Analitical(lambda2,z,5);
    uanali5 = [uanali5 abs(a2).^2];
    a3 = Analitical(lambda2,z,7.5);
    uanali75 = [uanali75 abs(a3).^2];
end

figure(Color='w');
subplot(1, 3, 1); plot(z2,uanali1); xlabel("z[m]"); ylabel("I analytical [a.u.]"); title("a = 1um");
subplot(1, 3, 2); plot(z2,uanali5); xlabel("z[m]"); ylabel("I analytical [a.u.]"); title("a = 5um");
subplot(1, 3, 3); plot(z2,uanali75); xlabel("z[m]"); ylabel("I analytical [a.u.]"); title("a = 7.5um");





