function [nrings,ringnames] = ReadRingsAndNames(infofile,hs_range)

nhs = length(hs_range);
T = readtable(infofile);
c = 0;
for i = 1:size(T,1)
    if T.heightstep(i) == hs_range(1)
        c = c+1;
    end
end
nrings = c;
for io = 1:nhs
    hs = hs_range(io);
    loopind = find(T.heightstep==hs);
    for c = 1:nrings
        scandir = T.scandir{loopind(c)};
        ringNo = T.ring(loopind(c));
        this_suffix = T.suffix(loopind(c));
        if this_suffix == 0
            scansuffix = '';
        else
            scansuffix = ['_' num2str(this_suffix)];
        end
        ringnames{io,ringNo} = [scandir scansuffix];
%         if isnumeric(T.addstring)
%             addstring = '';
%         else
%             addstring = T.addstring{loopind(c)};
%         end
%         ringnames{io,c} = [scandir scansuffix addstring];
    end
end


end

