clear;
close all;
clc;
%% parameter
load('object_grayscale');
F1 =imresize(F1,[512,512]);[n,m]= size(F1);
E=sum(sum(F1));El=0.5*E;
lamda=532e-6;k=2*pi/lamda;
dh=0.00374;
F=abs(sqrt(F1));F=padarray(F,[n/4,m/4]);[nn,mm]=size(F);
phi=exp(1i*2*pi*rand(nn,mm));
amp=rand(nn,mm);
%% band-limitation
bandlim_spe=padarray(ones(nn/2,mm/2),[nn/4,mm/4]);
bandlim_in=ones(n,m);
bandlim_in=padarray(bandlim_in,[n/4,m/4]);
bandlim_ou=ones(nn,mm)-bandlim_in;
incident=bandlim_spe;
%% first iteration
loop=100;
RMSE=zeros(loop,1);
figure
for i=2:loop 
   amp=bandlim_in.*F+bandlim_ou.*amp;
   E1=amp.*phi;
   sigma=i/loop;
   E2=fftshift(fft2(fftshift(E1)));
   E2_ave1=sqrt((E+El)*incident.^2/sum(sum(incident.^2)));
   E2_ave2=sqrt((E+El)*abs(E2).^2/sum(sum(abs(E2).^2)));
   E2_k=E2_ave1.*exp(1i*angle(E2));
   E2_o=E2_ave2.*exp(1i*angle(E2));
   E2_soft= (1-sigma)*E2_o+sigma*E2_k;
   es=fftshift(ifft2(fftshift(E2_soft)));
   amp=abs(es);
   amp_in=bandlim_in.*amp; amp_ou=bandlim_ou.*amp;
   amp=sqrt(E*(amp_in.^2)/sum(sum(amp_in.^2)))+sqrt(El*(amp_ou.^2)/sum(sum(amp_ou.^2)));
   I=amp((nn/2-n/2)+1:(nn/2+n/2),(mm/2-m/2)+1:(mm/2+m/2)).^2;
   I=E*I/sum(sum(I));
   P=mod(angle(es),2*pi);
   P=P((nn/2-n/2)+1:(nn/2+n/2),(mm/2-m/2)+1:(mm/2+m/2));
   imshow(I);
   Diff=double(I)-double(F1);
   MSE=gather(sum(Diff(:).^2)/numel(I));
   RMSE(i,1)=sqrt(MSE);
   phi=exp(1i*angle(es));
   phi_in=angle(phi).*bandlim_in;
end
phase=angle(phi);
An=angle(E2_k);
hologram=incident.*exp(1i*An);
%% output
Rec=fftshift(ifft2(fftshift(hologram)));
I=abs(Rec).^2;
I=I((nn/2-n/2)+1:(nn/2+n/2),(mm/2-m/2)+1:(mm/2+m/2));
I=E*(I/sum(sum(I)));
figure,imshow(I);


