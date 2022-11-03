function [angles,ip180] = ReadAngles(filename,h5AnglePath)

find_index = @(v,val) find(abs(v-val) == min(abs(v-val)),1);

angles = abs(h5read(filename,h5AnglePath));
ip360 = find_index(angles,360+angles(1));% index of projection with angle 360 degrees
ip180 = floor(ip360/2); % index of projection with angle 180 degrees
angles = angles(1:ip180);
if mean(angles(:))>180
    angles = angles-180;
end

end
