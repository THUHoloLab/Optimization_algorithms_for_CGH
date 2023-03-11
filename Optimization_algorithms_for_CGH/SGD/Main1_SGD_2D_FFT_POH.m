close all;clear;clc;tic
%% parameter
slm.pix = 3.74e-3;
slm.Nx = 384; slm.Ny = 384; 
opt.Nx = 2*slm.Nx; opt.Ny = 2*slm.Ny;
dh = 3.74e-3;
%% illumination pattern at SLM
slm.window = zeros(opt.Nx,opt.Ny);
slm.window( (opt.Nx/2-slm.Nx/2)+1 : (opt.Nx/2+slm.Nx/2) , (opt.Ny/2-slm.Ny/2)+1 : (opt.Ny/2+slm.Ny/2)) = ones(slm.Nx, slm.Ny);
opt.source = slm.window;
opt.in = ones(opt.Nx*2/3,opt.Ny*2/3);
opt.in = padarray(opt.in,[opt.Nx/6,opt.Ny/6]);
opt.ou = 1-opt.in;
%% constructe object matrix
load('object_grayscale');
obj = imresize(F1, [opt.Nx*2/3, opt.Ny*2/3]);
obj = padarray(obj,[opt.Nx/6,opt.Ny/6]);
Masks = obj;
E=sum(Masks(:));
El=0.5*E; Es=1.5*E;
F1=Masks((opt.Nx/2-opt.Nx/3)+1:(opt.Nx/2+opt.Nx/3),(opt.Ny/2-opt.Ny/3)+1:(opt.Ny/2+opt.Ny/3));
%% generation of the starting value
[LX,LY]=size(Masks);
pha=2*pi*rand(LX,LY);
pha=exp(1i*pha);
amp=opt.in.*sqrt(Masks)+opt.ou.*rand(LX,LY);
hologram=ifftshift(ifft2(ifftshift(amp .* pha)));
pha=angle(hologram);
%% optimization
times=300;
RMSE=zeros(times,1);
figure
for i=2:times
    [loss, df] = Gradient_L2_FFT_POH(pha, opt.source, opt.in, opt.Nx, opt.Ny, Masks);
    [updates, state] = Optmization_SGD_ADAM(df, []);
    pha=pha-updates;
%% show result
objectField=opt.source.*exp(1i.*pha);
imagez=fftshift(fft2(fftshift(objectField)));  
amp=abs(imagez);
I=amp((opt.Nx/2-opt.Nx/3)+1:(opt.Nx/2+opt.Nx/3),(opt.Ny/2-opt.Ny/3)+1:(opt.Ny/2+opt.Ny/3)).^2;
I=E*I/sum(sum(I));
P=mod(angle(imagez),2*pi);
P=P((opt.Nx/2-opt.Nx/3)+1:(opt.Nx/2+opt.Nx/3),(opt.Ny/2-opt.Ny/3)+1:(opt.Ny/2+opt.Ny/3));
imshow(I);
Diff=double(I)-double(F1);
MSE=gather(sum(Diff(:).^2)/numel(I));
RMSE(i,1)=sqrt(MSE);
phi_in=angle(imagez).*opt.in;
end
rec.phase = reshape(pha, [opt.Nx, opt.Ny]);
hologram = exp(1i.*rec.phase);
%% reconstruction
Rec = fftshift(fft2(fftshift(opt.source.*hologram)));  
I = abs(Rec).^2;
I=I((opt.Nx/2-opt.Nx/3)+1:(opt.Nx/2+opt.Nx/3),(opt.Ny/2-opt.Ny/3)+1:(opt.Ny/2+opt.Ny/3));
I= E*I/sum(sum(I));
Phase = mod(angle(Rec),2*pi);
figure,imshow(I);
toc