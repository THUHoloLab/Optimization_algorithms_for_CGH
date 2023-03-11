function [loss, df ] = Gradient_L2_3D_AS_POH( phase, source, Nx, Ny,  mask, HTrans)
% This function calculate the loss function and the gradient with repect to
% the phase-only hologram. Angular spectrum theory is utilized for propagation.
slice=size(mask,3);
df = zeros(Nx, Ny);
loss = 0; 
V = mask;
mass2 = sum(V(:));
phase = reshape(phase, [Nx, Ny]);
%% compute loss
    objectField = source.*exp(1i*phase);
    for s=1:slice
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans(:,:,s));
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    I = abs(imagez).^2;
    mass1 = sum(I(:));
    I = mass2*I/mass1;
    diffh = (I-V(:,:,s)).^2;
    L2 = sum(diffh(:)); 
    loss = loss + L2; 
%% Compute gradient 
    temph = 2*(I-V(:,:,s));
    temph = mass1*temph/mass2;
    temph = 2*imagez.*temph;
    temph = fftshift(fft2(fftshift(temph))).*HTrans(:,:,s);
    temph = ifftshift(ifft2(ifftshift(temph)));
    df = df + temph;
    df=source.*df;
    df = - real(df).*sin(phase) + imag(df) .* cos(phase);
    end
end
