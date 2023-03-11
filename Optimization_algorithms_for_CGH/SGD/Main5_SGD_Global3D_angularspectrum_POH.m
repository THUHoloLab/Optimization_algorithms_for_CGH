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
load('object_3D');
% F1=F1(:,:,1:3);
F1 = imresize(F1, [slm.Nx, slm.Ny,[]]);
% F1=F1(:,:,1);
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
    depth=20+50*(s-1)/(slice-1);
    lim_m=((2*delta_m*depth).^2+1).^(-1/2)./opt.lambda;
    lim_n=((2*delta_n*depth).^2+1).^(-1/2)./opt.lambda;
    bandlim_m=(lim_m-abs(fx));
    bandlim_n=(lim_n-abs(fy));
    bandlim_m=imbinarize(bandlim_m,0);
    bandlim_n=imbinarize(bandlim_n,0);
    bandlim_AS=bandlim_m.*bandlim_n;
    HTrans(:,:,s) = bandlim_AS.*exp(1i*opt.k*sqrt(1-(opt.lambda*fy).^2-(opt.lambda*fx).^2)*depth);
end
% figure,imshow(F1(:,:,1));
%% generation of the starting value
hologram = zeros(opt.Nx,opt.Ny);
for s=1:slice
    pha=exp(1i*2*pi*rand(opt.Nx,opt.Ny));
    amp=sqrt(Masks(:,:,s));
    wavefront=fftshift(fft2(fftshift(amp .* pha))) .* HTrans(:,:,s);
    wavefront=ifftshift(ifft2(ifftshift(wavefront)));
    hologram=hologram+wavefront;
end
Pha=angle(opt.source.*hologram);
%% optimization
times=100;
RMSE=zeros(times,1);
I=zeros(slm.Nx, slm.Ny,slice);
figure(1)
for i=2:times
    [loss, df] = Gradient_L2_3D_AS_POH(Pha, opt.source, opt.Nx, opt.Ny, Masks, HTrans);
    [updates, state] = Optmization_SGD_ADAM(df, []);
    Pha=Pha-updates;
%% show result
    objectField=opt.source.*exp(1i.*Pha);
    for s=1:slice
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans(:,:,s));
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    amp=abs(imagez);
    I(:,:,s)=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    end
    I=E*I/sum(I(:));
    imshow(I(:,:,1));
    Diff=double(I)-double(F1);
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
end
rec.phase = reshape(Pha, [opt.Nx, opt.Ny]);
rec.phase = exp(1i.*rec.phase);
%% reconstruction
f = figure(2);
for s=1:2*slice
%     if mod(s,2)== 1
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
