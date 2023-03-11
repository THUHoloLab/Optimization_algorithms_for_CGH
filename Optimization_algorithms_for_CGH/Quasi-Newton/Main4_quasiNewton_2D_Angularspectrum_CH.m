close all;clear;clc;tic
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
%% optimization initialization
matlab_options = optimoptions('fmincon','GradObj','on', ...
    'algorithm','interior-point','Hessian','lbfgs','FunValCheck','on','MaxFunEvals', 500 ,'MaxIter', 10,...
    'TolX', 1e-20, 'TolFun', 1e-15);
lb = -inf(opt.Nx*opt.Ny, 1);
ub = inf(opt.Nx*opt.Ny, 1);
times=30;
RMSE=zeros(times,1);
figure
pha=2*pi*rand(opt.Nx,opt.Ny);
Masks = obj;
E=sum(Masks(:));
f = @(x)Gradient_L2_AS_CH(x, opt.source, opt.Nx, opt.Ny, Masks, HTrans);
for i=2:times
    [pha, loss] = fmincon(f,pha,[],[],[],[],lb,ub,[],matlab_options);
%% show reconstructions
    objectField = sqrt(Masks).*exp(1i.*pha);
    imagez = fftshift(fft2(fftshift(objectField))) .* HTrans;
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    imagez = opt.source.*imagez;
    Am=sqrt(E*(abs(imagez).^2)./sum(sum(abs(imagez).^2)));
    wavefront= Am.*exp(1i*angle(imagez));
    recz = fftshift(fft2(fftshift(imagez))) .* conj(HTrans);
    recz = ifftshift(ifft2(ifftshift(recz)));
    amp=abs(recz);
    I=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
    I = E*I/sum(sum(I));
    imshow(I);
    Diff=double(I)-double(F1);
    MSE=gather(sum(Diff(:).^2)/numel(I));
    RMSE(i,1)=sqrt(MSE);
end
hologram=wavefront;
%% reconstruction
f = figure(2);
Masks = obj;
objectField = opt.source.* hologram;
rec = fftshift(fft2(fftshift(objectField))) .* conj(HTrans);
rec = ifftshift(ifft2(ifftshift(rec)));
amp=abs(rec);
I=amp(opt.Nx/4+1:opt.Nx*3/4,opt.Ny/4+1:opt.Ny*3/4).^2;
I = E*I/sum(sum(I));
toc
imwrite(I,'QN_AS_3D_CH.bmp');
imshow(I);        
