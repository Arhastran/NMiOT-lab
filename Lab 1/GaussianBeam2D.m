%% Function provided for labs 

function [Uout] = GaussianBeam2D(A0,z,z0,Nx,dx,lambda)
% creates gaussian Beam

%global x2
x = (-Nx/2:Nx/2-1)*dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%this is for different sign conversion for Saleh and Goodman
z = -z;

mone = ones(1,max(size(x)));
% z0 - Rayleigh Range 
k = 2*pi/lambda;
W0 = sqrt(lambda*z0/pi);
%czynnik przesuniecia fazowy
ksiz= atan(z/z0);
%promien fazy
Rz = z*(1+(z0/z)^2);
%srednica
Wz = W0*sqrt(1+(z/z0)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% REM X2 AND SHIFTED TO GLOBAL
x2 = mone'*(x.^2)+(x'.^2)*mone;

if z == 0 % waist 
    Uout = A0*W0/Wz*exp(- x2/Wz^2  ).*exp( -1i*k*z  + 1i*ksiz );  
else
    %!!!!!!!!!!!!!
    % this is good equation
   Uout = A0*W0/Wz*exp(- x2/Wz^2  ).*exp( -1i*k*z - 1i*k*x2/(2*Rz) + 1i*ksiz );  
    % this is wrong equation
%     Uout = A0*exp(- x2/Wz^2  ).*exp( -1i*k*z - 1i*k*x2/(2*Rz) + 1i*ksiz ); 
end