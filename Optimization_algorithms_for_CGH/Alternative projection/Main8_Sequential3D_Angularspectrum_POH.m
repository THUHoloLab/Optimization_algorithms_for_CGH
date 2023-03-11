close all;clear;clc;
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
load('Object_3D_binary2');
F1 = imresize(F1, [slm.Nx, slm.Ny,[]])+0.001;
obj = padarray(F1,[slm.Nx/2,slm.Ny/2,[]]);
Masks = obj;
E=sum(Masks(:));
opt.source=sqrt(E*opt.source/sum(opt.source(:)));
figure,imshow(Masks(:,:,3));
%% convolution kernel (bandlimited )
slice=size(Masks,3);
HTrans=zeros(opt.Nx,opt.Ny,slice);
[fx,fy]=meshgrid(linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Ny),linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Nx));
Sm=opt.Nx*dh;Sn=opt.Ny*dh;
delta_m=(2*Sm).^(-1);delta_n=(2*Sn).^(-1);
for s=1:slice
    if s==1
        depth=10;
    else
        depth=10;
    end
    lim_m=((2*delta_m*depth).^2+1).^(-1/2)./opt.lambda;
    lim_n=((2*delta_n*depth).^2+1).^(-1/2)./opt.lambda;
    bandlim_m=(lim_m-abs(fx));
    bandlim_n=(lim_n-abs(fy));
    bandlim_m=imbinarize(bandlim_m,0);
    bandlim_n=imbinarize(bandlim_n,0);
    bandlim_AS=bandlim_m.*bandlim_n;
    HTrans(:,:,s) = bandlim_AS.*exp(1i*opt.k*sqrt(1-(opt.lambda*fy).^2-(opt.lambda*fx).^2)*depth);
end
%% optimization
times=30;
RMSE=zeros(times,1);
I=zeros(slm.Nx, slm.Ny, slice);
pha=exp(1i*2*pi*rand(opt.Nx,opt.Ny));
figure(1)
for i=2:times
    for s=1:slice
        mask=Masks(:,:,slice-s+1);
%         amp=sqrt(E.*mask/sum(mask(:)));
        amp=mask;
        wavefront=fftshift(fft2(fftshift(amp .* pha))) .* HTrans(:,:,slice-s+1);
        wavefront=ifftshift(ifft2(ifftshift(wavefront)));
        pha=exp(1i*angle(wavefront));
    end
    objectField=opt.source.*exp(1i.*pha);
    for s=1:slice
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans(:,:,s));
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    Amp=abs(imagez);
    Pha=angle(imagez);
    mask=Masks(:,:,s);
%     amp=sqrt(E.*mask/sum(mask(:)));
    amp=mask;
    objectField = amp.*exp(1i.*Pha);
    I(:,:,s)=Amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    end
    pha=Pha;
    I=E*I/sum(I(:));
    imshow(I(:,:,2));
    Diff=double(I)-double(F1);
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
end
hologram = reshape(pha, [opt.Nx, opt.Ny]);
hologram = exp(1i.*hologram);
%% reconstruction
f = figure(2);
for s=1:2*slice
    if s<=slice
        subplot(2,3,s)
        imshow(F1(:,:,s));
    else
        subplot(2,3,s)       
        axis off; axis image; colormap gray;     
        imshow(I(:,:,s-slice));
    end
end
for s=1:slice
imwrite(I(:,:,s),[num2str(s),'Global_3D_POH.bmp']);
end
