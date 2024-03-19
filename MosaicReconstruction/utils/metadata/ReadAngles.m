function [angles,ip180] = ReadAngles(filename,h5AnglePath)

find_index = @(v,val) find(abs(v-val) == min(abs(v-val)),1);

angles = h5read(filename,h5AnglePath);
angles = abs(angles - angles(1));
ip360 = find_index(angles,360);% index of projection with angle 360 degrees
ip180 = floor(ip360/2); % index of projection with angle 180 degrees
angles = angles(1:ip180);

end
