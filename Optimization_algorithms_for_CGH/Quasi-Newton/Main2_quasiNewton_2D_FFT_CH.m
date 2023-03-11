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
%% Generate Dummy Data
load('object_grayscale');
obj = imresize(F1, [opt.Nx, opt.Ny]);
Masks=obj; E=sum(Masks(:));
%% generation of the starting value
[LX,LY]=size(Masks);
pha=2*pi*rand(LX,LY);
%% optimization
matlab_options = optimoptions('fmincon','GradObj','on', ...
    'algorithm','interior-point','Hessian','lbfgs','FunValCheck','on','MaxFunEvals', 500 ,'MaxIter', 300,...
    'TolX', 1e-20, 'TolFun', 1e-15);
lb = -inf(opt.Nx*opt.Ny, 1);
ub = inf(opt.Nx*opt.Ny, 1);
f = @(x)Gradient_L2_FFT_CH(x, opt.source, opt.Nx, opt.Ny, Masks);
times=10;
RMSE=zeros(times,1);
figure
for i=2:times
[phase, loss] = fmincon(f,pha,[],[],[],[],lb,ub,[],matlab_options);
%% show result
objectField = sqrt(Masks).*exp(1i.*phase);
imagez = opt.source.*fftshift(fft2(fftshift(objectField)));  
imagez = ifftshift(ifft2(ifftshift(imagez)));  
I = abs(imagez).^2;
I = E*I/sum(sum(I));
P = mod(angle(imagez),2*pi);
imshow(I);
RMSE(i,1)=sqrt(loss/numel(Masks));
pha=phase;
end
rec.phase = reshape(phase, [opt.Nx, opt.Ny]);
objectField = sqrt(Masks).*exp(1i.*rec.phase);
hologram = opt.source.*fftshift(fft2(fftshift(objectField)));  
%% reconstruction
Rec = ifftshift(ifft2(ifftshift(hologram)));  
I=abs(Rec).^2;
I=E*I/sum(sum(I));
imwrite(I,'CH.bmp');
Phase=mod(angle(Rec),2*pi);
figure,imshow(I);
toc