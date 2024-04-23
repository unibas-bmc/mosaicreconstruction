function [tzpos] = ReadTzPos(parfile)
myStr1 = '#tomo1tz = '; lenStr1 = length(myStr1);

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
tzpos = str2double(newStr(1));
end

