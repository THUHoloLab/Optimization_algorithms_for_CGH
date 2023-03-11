function [loss, df] = Gradient_L2_AS_CH( phase, source, Nx, Ny,  mask, HTrans)
% 
df = zeros(Nx, Ny);
loss = 0; 
V = mask;
mass2 = sum(V(:));
% phase = reshape(phase, [Nx, Ny]);

objectField = sqrt(V).*exp(1i*phase);
    imagez = fftshift(fft2(fftshift(objectField))) .* HTrans;
    imagez = ifftshift(ifft2(ifftshift(imagez))); 
    imagez = source.* imagez;
    imagez = fftshift(fft2(fftshift(imagez))) .* conj(HTrans);
    imagez = ifftshift(ifft2(ifftshift(imagez)));    
    I = abs(imagez.^2);
    mass1 = sum(I(:));
    I = mass2*I/mass1;
    diffh = (I-V).^2;
    L2 = sum(sum(diffh)); 
    loss = loss + L2; 
%% Compute gradient 
    temph = 2*(I-V);
    temph = mass1*temph/mass2;
    temph = 2*imagez.*temph;
    temph = fftshift(fft2(fftshift(temph))).* HTrans;
    temph = ifftshift(ifft2(ifftshift(temph)));
    temph = temph.*source;
    temph = fftshift(fft2(fftshift(temph))).* conj(HTrans);
    temph = ifftshift(ifft2(ifftshift(temph)));
    df = sqrt(V).*temph;
    df = real(- real(df).*sin(phase) + imag(df) .* cos(phase));
end
