close all;clear;clc;
%% parameter
slm.pix = 3.74e-3;
slm.Nx = 512; slm.Ny = 512;
opt.Nx = 2*slm.Nx; opt.Ny = 2*slm.Ny;
dh = 3.74e-3;
opt.lambda = 532e-6;
opt.k=2*pi/opt.lambda;
%% illumination pattern at SLM
opt.source = 1;
slm.window = zeros(opt.Nx,opt.Ny);
slm.window( (opt.Nx/2-slm.Nx/2)+1 : (opt.Nx/2+slm.Nx/2) , (opt.Ny/2-slm.Ny/2)+1 : (opt.Ny/2+slm.Ny/2)) = ones(slm.Nx, slm.Ny);
opt.source = opt.source.*slm.window;
%% construct object matrix
load('Object_3D_grayscale');
F1 = imresize(F1, [slm.Nx, slm.Ny,[]]);
obj = padarray(F1,[slm.Nx/2,slm.Ny/2,[]]);
%% convolution kernel (bandlimited )
slice=size(obj,3);
HTrans=zeros(opt.Nx,opt.Ny,slice);
[fx,fy]=meshgrid(linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Ny),linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Nx));
Sm=opt.Nx*dh;Sn=opt.Ny*dh;
delta_m=(2*Sm).^(-1);delta_n=(2*Sn).^(-1);
for s=1:slice
    depth=20+2*(s-1)/(slice-1);
    lim_m=((2*delta_m*depth).^2+1).^(-1/2)./opt.lambda;
    lim_n=((2*delta_n*depth).^2+1).^(-1/2)./opt.lambda;
    bandlim_m=(lim_m-abs(fx));
    bandlim_n=(lim_n-abs(fy));
    bandlim_m=imbinarize(bandlim_m,0);
    bandlim_n=imbinarize(bandlim_n,0);
    bandlim_AS=bandlim_m.*bandlim_n;
    HTrans(:,:,s) = bandlim_AS.*exp(1i*opt.k*sqrt(1-(opt.lambda*fy).^2-(opt.lambda*fx).^2)*depth);
end
%% optimization initialization
times=100;
hologram=zeros(opt.Nx, opt.Ny);
RMSE=zeros(times,1);
figure
for s=1:slice 
    pha=2*pi*rand(opt.Nx,opt.Ny);
    Masks = obj(:,:,s);
    E=sum(Masks(:));
    hTrans = HTrans(:,:,s);
    for i=2:times
%% show reconstructions
    objectField = sqrt(Masks).*exp(1i.*pha);
    imagez = fftshift(fft2(fftshift(objectField))) .* hTrans;
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    imagez = opt.source.*imagez;
    Am=sqrt(E*(abs(imagez).^2)./sum(sum(abs(imagez).^2)));
    wavefront= Am.*exp(1i*angle(imagez));
    recz = fftshift(fft2(fftshift(wavefront))) .* conj(hTrans);
    recz = ifftshift(ifft2(ifftshift(recz)));
    pha= angle(recz);
    amp=abs(recz);
    I=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    I = E*I/sum(sum(I));
    imshow(I);
    Diff=double(I)-double(F1(:,:,s));
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
    end
    hologram=hologram+wavefront;
end
hologram=angle(hologram);
%% reconstruction
f = figure(2);
E=sum(obj(:));
for s=1:2*slice
%     if mod(s,2)== 1
    if s<=slice
        subplot(2,3,s)
        imshow(F1(:,:,s));
    else
        subplot(2,3,s)
        Masks = obj(:,:,s-slice);
        objectField = opt.source.* exp(1i.*hologram);
        rec = fftshift(fft2(fftshift(objectField))) .* conj(HTrans(:,:,s-slice));
        rec = ifftshift(ifft2(ifftshift(rec)));
        amp=abs(rec);
        I=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
        I = E*I/sum(sum(I));
        imwrite(I,[num2str(s),'super_3D_POH.bmp']);
        imshow(I);        
    end
end
