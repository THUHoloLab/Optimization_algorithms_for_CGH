close all;clear;clc;tic
%% parameter
slm.pix = 3.74e-3;
slm.Nx = 256; slm.Ny = 256;
opt.Nx = 2*slm.Nx; opt.Ny = 2*slm.Ny;
dh = 3.74e-3;
%% illumination pattern at SLM
opt.source = 1;
slm.window = zeros(opt.Nx,opt.Ny);
slm.window( (opt.Nx/2-slm.Nx/2)+1 : (opt.Nx/2+slm.Nx/2) , (opt.Ny/2-slm.Ny/2)+1 : (opt.Ny/2+slm.Ny/2)) = ones(slm.Nx, slm.Ny);
opt.source = opt.source.*slm.window;
%% construct object matrix
load('object_grayscale');
obj = imresize(F1, [opt.Nx, opt.Ny]);
Masks=obj; E=sum(Masks(:));
%% generation of the starting value
[LX,LY]=size(Masks);
pha=2*pi*rand(LX,LY);
%% optimization initialization
times=300;
RMSE=zeros(times,1);
figure
for i=2:times
[loss, df] = Gradient_L2_FFT_CH( pha, opt.source, opt.Nx, opt.Ny, Masks);
[updates, state] = Optmization_SGD_ADAM(df, []);
pha=pha-updates;
% objectField = sqrt(Masks).*exp(1i.*pha);
% imagez = opt.source.*fftshift(fft2(fftshift(objectField)));  
% imagez = ifftshift(ifft2(ifftshift(imagez)));  
% I = abs(imagez).^2;
% I = E*I/sum(sum(I));
% P = mod(angle(imagez),2*pi);
% imshow(I);
% RMSE(i,1)=sqrt(loss/numel(Masks));
end
rec.phase = reshape(pha, [opt.Nx, opt.Ny]);
objectField = sqrt(Masks).*exp(1i.*rec.phase);
hologram = opt.source.*fftshift(fft2(fftshift(objectField)));
toc
%% reconstruction
Rec = ifftshift(ifft2(ifftshift(hologram)));  
I=abs(Rec).^2;
I=E*I/sum(sum(I));
imwrite(I,'Global_3D_POH.bmp');
Phase=mod(angle(Rec),2*pi);
figure,imshow(I);
toc