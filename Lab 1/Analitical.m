function [out] = Analitical(lambda,z,a,dx)

[Ny,Nx] = size(z);
dfx = 1/Nx/dx;
dfy = 1/Ny/dx;
fx=(-Nx/2:Nx/2-1)*dfx;
fy=(-Ny/2:Ny/2-1)*dfy;
[FX,FY] = meshgrid(fx,fy);

k = 2*pi/lambda;
Func = @(z) z.*(exp(i*k*z)./z-((exp(i.*k.*sqrt(z.^2+a.^2)))/(sqrt(z.^2+a^2))));
out = Func(z);
end