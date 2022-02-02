
function len = SegmentLength(iHole)
% Copyright (C) 2016 - Universtiy of Montpellier
%
% About your license if you have any
%
% Version 0.0.1 - Denis Mottet - 1 August 2016

if size(iHole, 1) > size(iHole, 2)
    iHole = iHole';
end
len = diff([0,   find(diff(iHole) > 1 ),  length(iHole)]);

    
end
