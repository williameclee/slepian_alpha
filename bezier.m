% XYfine = BEZIER(XY)
%
% An attempt at smoothing a coastline by B-spline fitting.
%
% INPUT:
%
% XY        The set of points, make sure XY(end,:)=XY(1,:)
% N         Number of times this needs to be smoothed
%
% OUTPUT:
%
% XYfine    The smoothed curve
%
% http://www.me.cmu.edu/faculty1/shimada/gm98/project/ivan/project/
% Fujio Yamaguci "Curves and Surfaces in Computer Aided Geometric
% Design", Springer-Verlag, Berlin, 1988
% http://mathworld.wolfram.com/CubicSpline.html
%
% Last modified by fjsimons-at-alum.mit.edu, June 4th, 2004
% Last modified by charig-at-princeton.edu, April 24th, 2015
% Last modified by williameclee-at-arizona.edu, June 10th, 2024

function XYfine = bezier(XY, N)

    % Split the curve into cell of segments
    XYsplit = splitxy(XY);

    % Smooth each segment
    XYsplit = cellfun(@(x) splinesmoothing(x, N), ...
        XYsplit, 'UniformOutput', false);

    % Add NaNs between segments
    XYsplit = cellfun(@(x) [x; NaN NaN], ...
        XYsplit, 'UniformOutput', false);

    % Join the segments back together
    XYfine = cell2mat(XYsplit);
    XYfine = XYfine(1:end - 1, :);

    % end

end

function XYfine = splinesmoothing(XY, N)
    % Seperated the spline smoothing from the bezier function
    % To avoid bezier calling itself recursively

    % Remove duplicate vertices
    XY = removeduplicatevertices(XY);

    % Number of points of the input curve
    n = size(XY, 1);

    % Calculate cumulative geodesic distance between these points
    % This is the "knot vector"
    cumDist = cumsum( ...
        [0; grcdist(XY(1:end - 1, :), XY(2:end, :))]);
    cumDist = cumDist / cumDist(end);
    % New, finer, linear distance between the points
    cumDistFine = linspace(0, 1, n * N);

    % Smoothed result
    XYfine = spline(cumDist, XY', cumDistFine)';
end
