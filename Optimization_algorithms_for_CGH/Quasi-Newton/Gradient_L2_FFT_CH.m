function [loss, df] = Gradient_L2_FFT_CH( phase, source, Nx, Ny,  mask)

df = zeros(Nx, Ny);
loss = 0; 
V = mask;
mass2 = sum(V(:));
phase = reshape(phase, [Nx, Ny]);
% amplitude=abs(complex);
% phase=angle(complex);

objectField = sqrt(V).*exp(1i*phase);
    imagez = fftshift(fft2(fftshift(objectField))); 
    imagez = ifftshift(ifft2(ifftshift(source.*imagez)));
    I = abs(imagez.^2);
    mass1 = sum(I(:));
    I = mass2*I/mass1;
%     I = I/mass;
%     V = mask/sum(mask(:)); 
    diffh = (I-V).^2;
    L2 = sum(sum(diffh)); 
       %Total variation
        %Compute losses
    loss = loss + L2; 
        %Compute gradient 
    temph = 2*(I-V);
    temph = mass1*temph/mass2;
    temph = 2*imagez.*temph;
    temph = fftshift(fft2(fftshift(temph))); 
    temph = ifftshift(ifft2(ifftshift(source.*temph)));
    df = sqrt(V).*temph;
    df = real(- real(df).*sin(phase) + imag(df) .* cos(phase));
%     dfphase(:,1)=df(:);
    
loss = gather(real(loss));
end
