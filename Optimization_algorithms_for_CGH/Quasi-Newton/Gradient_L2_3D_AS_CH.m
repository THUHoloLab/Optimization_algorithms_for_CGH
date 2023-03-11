function [loss, df ] = Gradient_L2_3D_AS_CH( phase, source, Nx, Ny,  mask, HTrans)
% This function calculate the loss function and the gradient with repect to
% the phase-only hologram. Angular spectrum theory is utilized for propagation.
slice=size(mask,3);
df = zeros(Nx, Ny, slice);
recz=zeros(Nx, Ny, slice);
loss = 0; 
V = mask;
mass2 = sum(V(:));
wavefront=zeros(Nx, Ny);
% phase = reshape(phase, [Nx, Ny,slice]);
%% compute loss
for s=1:slice
    objectField = sqrt(V(:,:,s)).*exp(1i*phase(:,:,s));
    imagez = fftshift(fft2(fftshift(objectField))) .* HTrans(:,:,s);
    imagez = ifftshift(ifft2(ifftshift(imagez)));
    wavefront=wavefront+imagez;
end
    hologram = source.* wavefront;
for s=1:slice
    rec = fftshift(fft2(fftshift(hologram))) .* conj(HTrans(:,:,s));
    recz(:,:,s) = ifftshift(ifft2(ifftshift(rec)));  
end
    I = abs(recz).^2;
    mass1 = sum(I(:));
    I = mass2*I/mass1;
    diffh = (I-V).^2;
    L2 = sum(diffh(:)); 
    loss = L2; 
%% Compute gradient 
    temph = 2*(I-V);
    temph = mass1*temph/mass2;
    temph = 2*imagez.*temph;
    temph = fftshift(fft2(fftshift(temph))).*HTrans;
    temph = ifftshift(ifft2(ifftshift(temph)));
    temph = temph.*source;
    temph = fftshift(fft2(fftshift(temph))).* conj(HTrans);
    temph = ifftshift(ifft2(ifftshift(temph)));
    df=  sqrt(V).*temph;
    df = - real(df).*sin(phase) + imag(df) .* cos(phase);
end
