function [loss, df ] = Gradient_L2_FFT_POH( phase, source, field, Nx, Ny,  mask)

df = zeros(Nx, Ny);
loss = 0; 
V = mask;
mass2 = sum(V(:));
phase = reshape(phase, [Nx, Ny]);
% amplitude=abs(complex);
% phase=angle(complex);

objectField = source.*exp(1i*phase);
    imagez = fftshift(fft2(fftshift(objectField))); 
    imagez = field.*imagez;
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
    temph = temph.*field;
    temph = ifftshift(ifft2(ifftshift(temph)));
    df = df + temph;
    df=source.*df;
    dfphase = - real(df).*sin(phase) + imag(df) .* cos(phase);
    df = real(dfphase);
loss = gather(real(loss));
end
