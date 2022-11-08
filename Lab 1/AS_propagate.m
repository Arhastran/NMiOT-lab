%% Function provided for labs 

function uout = AS_propagate(uin, z, lambda, n0, dx)

[Ny,Nx] = size(uin);
dfx = 1/Nx/dx;
dfy = 1/Ny/dx;
fx=(-Nx/2:Nx/2-1)*dfx;
fy=(-Ny/2:Ny/2-1)*dfy;
[FX,FY] = meshgrid(fx,fy);

%FT of the input field
FTu = fftshift(fft2(uin));

% %generation of the transfer function
FZ = sqrt((n0/lambda).^2-FX.^2-FY.^2);
FZ(~isreal(FZ))=0;
TF = exp(1i*2*pi*z*FZ);

% multiplication with the transfer function
FTu = FTu.*TF;

% inverse FT
uout = ifft2(fftshift(FTu));

