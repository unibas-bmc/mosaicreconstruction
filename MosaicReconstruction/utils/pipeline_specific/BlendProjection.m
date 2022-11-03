function [im_blend] = BlendProjection(im,mask)
%BLENDPROJ Blends projection based on provided mask
%   Detailed explanation goes here

%im_blend = nanmean(im.*(1+0./mask),3);
im_blend = sum(im.*mask,3)./sum(mask,3);

if any(isnan(im_blend),'all')
    tmp = sum(im,3);
    im_blend(isnan(im_blend)) = tmp(isnan(im_blend));
end

end

