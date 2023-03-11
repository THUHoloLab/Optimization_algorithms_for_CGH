function [loss, df ] = Gradient_L2_AS_POH( phase, source, Nx, Ny,  mask, HTrans)
% This function calculate the loss function and the gradient with repect to
% the phase-only hologram. Angular spectrum theory is utilized for propagation.

df = zeros(Nx, Ny);
loss = 0; 
V = mask;
mass2 = sum(V(:));
phase = reshape(phase, [Nx, Ny]);
%% compute loss
    objectField = source.*exp(1i*phase);
    imagez = fftshift(fft2(fftshift(objectField))) .* conj(HTrans);
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    I = abs(imagez).^2;
    mass1 = sum(I(:));
    I = mass2*I/mass1;
    diffh = (I-V).^2;
    L2 = sum(sum(diffh)); 
    loss = loss + L2; 
%% Compute gradient 
    temph = 2*(I-V);
    temph = mass1*temph/mass2;
    temph = 2*imagez.*temph;
%     temph = temph.*field;
    temph = fftshift(fft2(fftshift(temph))).*HTrans;
    temph = ifftshift(ifft2(ifftshift(temph)));
    df = df + temph;
    df=source.*df;
    dfphase = - real(df).*sin(phase) + imag(df) .* cos(phase);
    df = real(dfphase);
end
