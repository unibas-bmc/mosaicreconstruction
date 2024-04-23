function shifted_image = subpixelshift(img,rowShift,colShift)
    img = fft2(img);
    [nr,nc]  = size(img);
    Nr = ifftshift(-round(nr/2):(floor(nr/2)-1));
    Nc = ifftshift(-round(nc/2):(floor(nc/2)-1));
    img = img.* exp(2i*pi*((-rowShift*Nr/nr).'-(colShift*Nc/nc)));
    %Greg = Greg * exp(1im*diffphase)
    %shifted_image = abs(ifft2(Greg));
    shifted_image = real(ifft2(img));
    %shifted_image = abs(shifted_image);
end 
