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
Nx = 2048; %size of matrix 
Ny = 2048; 
dx = lambda2/2; %
dy = lambda2/2; %
d = 40; %size of aperture (side or 2*radius)
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

% XD = Analitical(lambda2,uin1,300,dx);
% imagesc(x,y,abs(XD));


%Square aperture+signal 

% d2 = 24;
% d3 = 2024;
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
uout_theory = AS_propagate(gauss,z1+prop_dist,lambda2,n0,dx);
% for zero padding method 
new_theory = zeros(size(uout_theory)+2);
new_theory(2:end-1,2:end-1)=uout_theory;
%uout_normal = AS_propagate(gauss, prop_dist, lambda2, n0, dx);
diff = (uout_theory)-(uout1);
figure('Color','w');
subplot(1,2,1); plot(x,abs(gauss(:,Nx/2+1)).^2);title("Gaussian beam before propagation");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(1,2,2); plot(x,abs(uout_theory(:,Nx/2+1)).^2);title("Gaussian beam after propagation");xlabel('x [um]');ylabel('Intensity[a.u]');

figure('Color','w', 'Name', 'Comparision of intensity distribution and error bars');
subplot(2,2,1), plot(x,abs(uout_theory(:,Nx/2+1)).^2);title("Gaussian beam numerical");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,2,2), plot(x,abs(uout1(:,Nx/2+1)).^2);title("Intensity after diffraction and propagation");xlabel('x [um]');ylabel('Intensity[a.u]');
subplot(2,2,3), imagesc(x,y,abs(uout_theory)); title("Amp out numerical Gauss");xlabel('x [um]');ylabel('y [um]');colorbar;
subplot(2,2,4), imagesc(x,y,abs(uout1)); title("Amp out diffraction");xlabel('x [um]');ylabel('y [um]');colorbar;

newin1 = zeros(size(uin1)+2);
newin1(2:end-1,2:end-1)=uin1;

newin2 = zeros(size(uin2)+2);
newin2(2:end-1,2:end-1)=uin2;

uout1_zero = AS_propagate(newin1, z1+prop_dist, lambda2, n0, dx);
uout2_zero = AS_propagate(newin2, z1+prop_dist, lambda2, n0, dx);

newx = [0 x 0];
newy = [0 y 0];


 figure('Color','w', 'Name', 'Comparision of intensity distribution with and without zero padding');
 subplot(2,2,1); imagesc(x,y,abs(uout1_zero)); title("Zero Padding");xlabel('x [um]');ylabel('Intensity[a.u]');
 subplot(2,2,2); imagesc(x,y,abs(uout1)); title("Without 0 padding");xlabel('x [um]');ylabel('Intensity[a.u]');
 subplot(2,2,3); imagesc(x,y,abs(uout2_zero)); title("Zero Padding");xlabel('x [um]');ylabel('Intensity[a.u]');
 subplot(2,2,4); imagesc(x,y,abs(uout2)); title("Without 0 padding");xlabel('x [um]');ylabel('Intensity[a.u]');
% subplot(2,3,1), plot(x,abs(new_theory(:,Nx/2+1)).^2);title("Gaussian beam theory");xlabel('x [um]');ylabel('Intensity[a.u]');
% subplot(2,3,2), plot(x,abs(newout1(:,Nx/2+1)).^2);title("Gaussian beam for numerical propagation");xlabel('x [um]');ylabel('Intensity[a.u]');
% subplot(2,3,3), plot(x,abs(diff(:,Nx/2+1)).^2);title("Error difference");xlabel('x [um]');ylabel('Error[a.u.]');
% subplot(2,3,4), imagesc(newx,newy,abs(newout1)); title("Amp out theory");xlabel('x [um]');ylabel('y [um]');colorbar;
% subplot(2,3,5), imagesc(x,y,abs(uout1)); title("Amp out intensity");xlabel('x [um]');ylabel('y [um]');colorbar;
% subplot(2,3,6), imagesc(newx,newy,abs(new_theory-newout1)); title("Amp diff");xlabel('x [um]');ylabel('y [um]');colorbar;
% 








