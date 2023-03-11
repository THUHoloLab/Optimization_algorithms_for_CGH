clear;
close all;
clc;
tic
%% parameter
load('object_grayscale');
F1 =imresize(F1,[512,512]);[n,m]= size(F1);
E=sum(sum(F1));
lamda=532e-6;k=2*pi/lamda;
dh=0.00374;
F=abs(sqrt(F1));[nn,mm]=size(F);
phi=exp(1i*2*pi*rand(nn,mm));
amp=rand(nn,mm);
%% band-limitation
bandlim_spe=padarray(ones(nn/2,mm/2),[nn/4,mm/4]);
incident=bandlim_spe;
incident=sqrt(E*incident.^2/sum(sum(incident.^2)));
source=ones(nn,mm);
%% first iteration
loop=300;
RMSE=zeros(loop,1);
figure
for i=2:loop 
   amp=F;
   E1=amp.*phi;
   E2=fftshift(fft2(fftshift(E1)));
   E2_ave1=sqrt(E*incident.^2/sum(sum(E*incident.^2)));
   E2_k=E2_ave1.*exp(1i*angle(E2));
   es=fftshift(ifft2(fftshift(E2_k)));
  %% show reconstruction
%    E2_s=incident.*exp(1i*angle(E2));
%    eo=fftshift(ifft2(fftshift(E2_s)));
   amp=abs(es);
   amp=sqrt(E*(amp.^2)/sum(sum(amp.^2)));
   I=amp.^2;
   I=E*I/sum(sum(I));
%    imwrite(I,'POH1.bmp');
   imshow(I);
   Diff=double(I)-double(F1);
   MSE=gather(sum(Diff(:).^2)/numel(I));
   RMSE(i,1)=sqrt(MSE);
   phi=exp(1i*angle(es));
end
phase=angle(phi);
An=angle(E2_k);
hologram=incident.*exp(1i*An);
%% output
Rec=fftshift(ifft2(fftshift(hologram)));
I=abs(Rec).^2;
I=E*(I/sum(sum(I)));
Diff=double(I)-double(F1);
MSE=gather(sum(Diff(:).^2)/numel(I));
RMSE2=sqrt(MSE);
% imwrite(I,'POH2.bmp');
figure,imshow(I);
toc
