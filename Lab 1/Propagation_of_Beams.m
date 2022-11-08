%% Adam Ignaciuk 
% script for 1st lab in NMiOT laboratory (using scripts and functions from supervisors)
clear
close all

%constants 
lambda = 0.6*1e-9;
lambda2 = 0.5;

z1=10; % distance from the waist

%constants required for Gaussian Beam 
A0 = 1; %Amplitude 
%z0 = 3000; %Distance 
Nx = 2048; %?
Ny = 2048; 
dx = lambda2/2; %
dy = lambda2/2; %
d = 100;
z0=d/2; % Rayleigh range
prop_dist = 2000;
n0 = 1;


%Circle 
x=(-Nx/2+1:Nx/2)*dx; %x points 
y=(-Ny/2+1:Ny/2)*dy; %y points 
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
aperture = zeros(Ny,Nx);
aperture(R<d/2)=1;
uin1 = aperture;  %signal entering the setup 

%Square 

d2 = 784; %first edge 
d3 = 1264; %second edge 
aperturesquare = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if i>d2 && i<d3 && j>d2 && j<d3
        aperturesquare(i,j) = 1;
        else 
        aperturesquare(i,j) = 0;
        end

    end

end


uin2 = aperturesquare;



uout1 = AS_propagate(uin1, z1+prop_dist, lambda2, n0, dx);
uout2 = AS_propagate(uin2, z1+prop_dist, lambda2, n0, dx);


%zero padding
newout1 = zeros(size(uout1)+2);
newout1(2:end-1,2:end-1)=uout1;

newout2 = zeros(size(uout2)+2);
newout2(2:end-1,2:end-1)=uout2;

newx = [0 x 0];
newy = [0 y 0];

figure('Color', 'w','Name', 'Circular aperture and beam propagation after it');
subplot(2,2,1), imagesc(x,y,abs(uin1));title('amp at the input plane of circle');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,2),imagesc(x,y,angle(uin1));title('phase at the input plane of circle');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,3),imagesc(x,y,abs(uout1));title('amp at the output plane of circle');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,4), imagesc(x,y,angle(uout1));title('phase at the output plane of circle');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
figure('Color', 'w','Name', 'Rectangular aperture and beam propagation after it');
subplot(2,2,1), imagesc(x,y,abs(uin2));title('amp at the input plane of square');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,2),imagesc(x,y,angle(uin2));title('phase at the input plane of square');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,3), imagesc(x,y,abs(uout2));title('amp at the output plane of sqare');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;
subplot(2,2,4), imagesc(x,y,angle(uout2));title('phase at the output plane of square');xlabel('x [um]');ylabel('y [um]');axis image;colorbar;

gauss = GaussianBeam2D(A0,z1,z0,Nx,dx,lambda2);
uout_theory = GaussianBeam2D(A0,z1+prop_dist,z0,Nx,dx,lambda2);
new_theory = zeros(size(uout_theory)+2);
new_theory(2:end-1,2:end-1)=uout_theory;
%uout_normal = AS_propagate(gauss, prop_dist, lambda2, n0, dx);
diff = (new_theory)-(newout1);

figure('Color','w', 'Name', 'Comparision of intensity distribution and error bars');
subplot(2,3,1), plot(x,abs(uout_theory(:,Nx/2+1)).^2);title("Gaussian beam theory");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,3,2), plot(x,abs(uout1(:,Nx/2+1)).^2);title("Gaussian beam for numerical propagation");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,3,3), plot(newx,abs(diff(:,Nx/2+1)).^2);title("Error difference");xlabel('x [um]');ylabel('Error[a.u.]');
subplot(2,3,4), imagesc(x,y,abs(uout_theory)); title("Amp out theory");xlabel('x [um]');ylabel('y [um]');colorbar;
subplot(2,3,5), imagesc(x,y,abs(uout1)); title("Amp out intensity");xlabel('x [um]');ylabel('y [um]');colorbar;
subplot(2,3,6), imagesc(x,y,abs(new_theory-newout1)); title("Amp diff");xlabel('x [um]');ylabel('y [um]');colorbar;

figure('Color','w', 'Name', 'Comparision of intensity distribution and error bars[zero padding]');
subplot(2,3,1), plot(newx,abs(new_theory(:,Nx/2+1)).^2);title("Gaussian beam theory");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,3,2), plot(newx,abs(newout1(:,Nx/2+1)).^2);title("Gaussian beam for numerical propagation");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,3,3), plot(newx,abs(diff(:,Nx/2+1)).^2);title("Error difference");xlabel('x [um]');ylabel('Error[a.u.]');
subplot(2,3,4), imagesc(newx,newy,abs(newout1)); title("Amp out theory");xlabel('x [um]');ylabel('y [um]');colorbar;
subplot(2,3,5), imagesc(x,y,abs(uout1)); title("Amp out intensity");xlabel('x [um]');ylabel('y [um]');colorbar;
subplot(2,3,6), imagesc(newx,newy,abs(new_theory-newout1)); title("Amp diff");xlabel('x [um]');ylabel('y [um]');colorbar;




