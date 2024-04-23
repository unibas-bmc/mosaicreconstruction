function [fprojstack] = paganin_filter_stack(projstack,deltabeta,pixsize_mm,lambda_mm,det_dist_mm,splitflag)
if nargin<6
    splitflag = '';
end

if strcmp(splitflag,'split')
    tmp = factor(size(projstack,3));
    nb = tmp(1); 
    bs = size(projstack,3)/nb;
    fprintf(['Splitting stack into ' num2str(nb) ' blocks...\n'])
    
    fprojstack = zeros(size(projstack));
    for i = 1:nb
        this_zcrop = (i-1)*bs+1:i*bs;
        
        ft_projstack = fftshift(fftshift(fft2(projstack(:,:,this_zcrop)),1),2);

        [u,v] = ndgrid(linspace(-size(ft_projstack,1)/2,size(ft_projstack,1)/2-1,size(ft_projstack,1)),...
            linspace(-size(ft_projstack,2)/2,size(ft_projstack,2)/2-1,size(ft_projstack,2)));
        u = 2*pi*u/(pixsize_mm*size(ft_projstack,1));
        v = 2*pi*v/(pixsize_mm*size(ft_projstack,2));

        ft_projstack = ft_projstack./(1 + (deltabeta*det_dist_mm*lambda_mm/4/pi)*(u.^2 + v.^2));
        fprojstack(:,:,this_zcrop) = abs(ifft2(ifftshift(ifftshift(ft_projstack,1),2)));
    end
    
else
ft_projstack = fftshift(fftshift(fft2(projstack),1),2);

[u,v] = ndgrid(linspace(-size(ft_projstack,1)/2,size(ft_projstack,1)/2-1,size(ft_projstack,1)),...
    linspace(-size(ft_projstack,2)/2,size(ft_projstack,2)/2-1,size(ft_projstack,2)));
u = 2*pi*u/(pixsize_mm*size(ft_projstack,1));
v = 2*pi*v/(pixsize_mm*size(ft_projstack,2));

ft_projstack = ft_projstack./(1 + (deltabeta*det_dist_mm*lambda_mm/4/pi)*(u.^2 + v.^2));
fprojstack = abs(ifft2(ifftshift(ifftshift(ft_projstack,1),2)));
%fprojstack =  (deltabeta/2) * log(abs(fprojstack));
end
end

