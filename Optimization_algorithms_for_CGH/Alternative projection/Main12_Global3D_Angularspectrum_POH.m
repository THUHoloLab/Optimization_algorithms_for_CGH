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
opt.in = slm.window;
%% construct object matrix
load('Object_3D_binary');
F1 = imresize(F1, [slm.Nx, slm.Ny,[]]);
obj = padarray(F1,[slm.Nx/2,slm.Ny/2,[]]);
Masks = obj;
E=sum(Masks(:));
%% convolution kernel (bandlimited )
slice=size(Masks,3);
HTrans=zeros(opt.Nx,opt.Ny,slice);
[fx,fy]=meshgrid(linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Ny),linspace(-1/(2*slm.pix),1/(2*slm.pix),opt.Nx));
Sm=opt.Nx*dh;Sn=opt.Ny*dh;
delta_m=(2*Sm).^(-1);delta_n=(2*Sn).^(-1);
for s=1:slice
    depth=30+50*(s-1)/(slice-1);
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
times=100;
RMSE=zeros(times,1);
I=zeros(slm.Nx, slm.Ny, slice);
hologram = zeros(opt.Nx,opt.Ny);
pha=exp(1i*2*pi*rand(opt.Nx,opt.Ny,slice));
figure(1)
for i=2:times
    for s=1:slice
        amp=sqrt(Masks(:,:,s));
        wavefront=fftshift(fft2(fftshift(amp .* pha(:,:,s)))) .* HTrans(:,:,s);
        wavefront=ifftshift(ifft2(ifftshift(wavefront)));
        hologram=hologram+wavefront;
    end
    Pha=angle(hologram);
    objectField=opt.source.*exp(1i.*Pha);
    for s=1:slice
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans(:,:,s));
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    amp=abs(imagez);
    pha(:,:,s)=angle(imagez);
    I(:,:,s)=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    end
    I=E*I/sum(I(:));
    imshow(I(:,:,1));
    Diff=double(I)-double(F1);
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
end
hologram = reshape(Pha, [opt.Nx, opt.Ny]);
hologram = exp(1i.*hologram);
%% reconstruction
f = figure(2);
for s=1:2*slice
    if s<=slice
        subplot(2,5,s)
        imshow(F1(:,:,s));
    else
        subplot(2,5,s)       
        axis off; axis image; colormap gray;     
        imshow(I(:,:,s-slice));
    end
end
for s=1:slice
    imwrite(I(:,:,s),[num2str(s),'Global_3D_POH.bmp']);
end
% saveas(f,'Grlobal_3D_POH.bmp');