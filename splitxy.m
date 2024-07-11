%%  SPLITXY
%   Split NaN-separated segments of a curve into separate cells
%   Last modified by williameclee-at-arizona.edu, June 10th, 2024

function XYcell = splitxy(XY)
    [Xcell, Ycell] = polysplit(XY(:, 1), XY(:, 2));
    XYcell = cellfun(@(x, y) [x, y], Xcell, Ycell, 'UniformOutput', false);

end
