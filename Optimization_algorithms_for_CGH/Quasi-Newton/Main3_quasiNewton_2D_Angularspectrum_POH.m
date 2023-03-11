close all;clear;clc;tic
%% parameter
slm.pix = 3.74e-3;
slm.Nx = 512; slm.Ny = 512; 
opt.Nx = 2*slm.Nx; opt.Ny = 2*slm.Ny;
dh = 3.74e-3;
opt.lambda = 532e-6;
opt.k=2*pi/opt.lambda;
%% illumination pattern at SLM
slm.window = ones(slm.Nx, slm.Ny);
slm.window = padarray(slm.window,[slm.Nx/2,slm.Ny/2]);
opt.source = slm.window;
%% construct object matrix
load('object_grayscale');
F1 = imresize(F1, [slm.Nx, slm.Ny]);
obj = padarray(F1,[slm.Nx/2,slm.Ny/2]);
Masks = obj;
E=sum(Masks(:));
%% convolution kernel (band-limited angular spectrum)
depth = 25;
[fx,fy]=meshgrid(linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Ny),linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Nx));
Sm=opt.Nx*dh;Sn=opt.Ny*dh;
delta_m=(2*Sm).^(-1);delta_n=(2*Sn).^(-1);
lim_m=((2*delta_m*depth).^2+1).^(-1/2)./opt.lambda;
lim_n=((2*delta_n*depth).^2+1).^(-1/2)./opt.lambda;
bandlim_m=(lim_m-abs(fx));
bandlim_n=(lim_n-abs(fy));
bandlim_m=imbinarize(bandlim_m,0);
bandlim_n=imbinarize(bandlim_n,0);
bandlim_AS=bandlim_m.*bandlim_n;
HTrans = bandlim_AS.*exp(1i*opt.k*sqrt(1-(opt.lambda*fy).^2-(opt.lambda*fx).^2)*depth);
% figure,imshow(F1);
%% generation of the starting value
[LX,LY]=size(Masks);
pha=2*pi*rand(LX,LY);
pha=exp(1i*pha);
amp=sqrt(Masks);
hologram=fftshift(fft2(fftshift(amp .* pha))) .* HTrans;
hologram = ifftshift(ifft2(ifftshift(hologram)));
Pha=angle(hologram);
%% optimization
matlab_options = optimoptions('fmincon','GradObj','on', ...
    'algorithm','interior-point','Hessian','lbfgs','FunValCheck','on','MaxFunEvals', 500 ,'MaxIter', 10,...
    'TolX', 1e-20, 'TolFun', 1e-15);
lb=-inf(opt.Nx*opt.Ny, 1);
ub=inf(opt.Nx*opt.Ny, 1);
f=@(x)Gradient_L2_AS_POH(x, opt.source, opt.Nx, opt.Ny, Masks, HTrans);
times=30;
RMSE=zeros(times,1);
figure
for i=2:times
    [Pha, loss] = fmincon(f,Pha,[],[],[],[],lb,ub,[],matlab_options);
%% show result
    objectField=opt.source.*exp(1i.*Pha);
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans);
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    amp=abs(imagez);
    I=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    I=E*I/sum(sum(I));
    imshow(I);
    Diff=double(I)-double(F1);
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
end
rec.phase = reshape(Pha, [opt.Nx, opt.Ny]);
rec.phase = exp(1i.*rec.phase);
%% reconstruction
Rec = fftshift(fft2(fftshift(opt.source.*rec.phase))).* conj(HTrans); 
Rec = ifftshift(ifft2(ifftshift(Rec))); 
I = abs(Rec).^2;
% I=I((opt.Nx/2-opt.Nx/3)+1:(opt.Nx/2+opt.Nx/3),(opt.Ny/2-opt.Ny/3)+1:(opt.Ny/2+opt.Ny/3));
I= E*I/sum(sum(I));
imwrite(I,'QN_2D_AS_POH.bmp');
Phase = mod(angle(Rec),2*pi);
toc
figure,imshow(I);