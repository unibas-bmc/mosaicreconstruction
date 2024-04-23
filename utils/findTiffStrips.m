function [strips, position] = findTiffStrips(yrange, stripwidth)
%FINDTIFFSTRIPS Range of strips in Tiff that contain given y-index range
%   For the given range of lines, and the given width of strips,
%   find the range of strips needed to cover this range of lines.
%   Also return the index range within these strips.
firstStrip = floor((min(yrange) - 1) / stripwidth) + 1;
lastStrip = floor((max(yrange) - 1) / stripwidth) + 1;
strips = firstStrip:lastStrip;
firstLine = min(yrange) - (firstStrip - 1) * stripwidth;
lastLine = max(yrange) - (firstStrip - 1) * stripwidth;
position = firstLine:lastLine;
end

