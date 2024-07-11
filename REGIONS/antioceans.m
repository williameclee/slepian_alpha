function varargout = antioceans(varargin)
    XY = oceans(varargin{:});
    
    % Get the inverted polygon
    ployXY = polyshape(XY(:, 1), XY(:, 2));
    polyWorld = polyshape([0, 360, 360, 0], [-90, -90, 90, 90]);
    polyXY = subtract(polyWorld, ployXY);
    XY = polyXY.Vertices;
    XY = closecoastline(XY);

    if nargout == 0
        figure
        plot(XY(:, 1), XY(:, 2), 'k-'); axis equal; grid on
    else
        varargout = {XY};
    end
end