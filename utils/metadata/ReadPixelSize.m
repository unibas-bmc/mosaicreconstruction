function [pixsize_mm] = ReadPixelSize(parfile)
myStr1 = 'IMAGE_PIXEL_SIZE_1 = '; lenStr1 = length(myStr1);

fid=fopen(parfile);
i = 1;
tline = fgetl(fid);
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    if contains(tline,myStr1)
        break
    end
end
fclose(fid);

tline = tline(lenStr1+1:end);
newStr = split(tline,' ');
pixsize_um = str2double(newStr(1));
pixsize_mm = pixsize_um*1e-3;
end

