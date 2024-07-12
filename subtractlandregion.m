function landRegion = ...
        subtractlandregion(landRegion, lonlim, latlim, varargin)
    p = inputParser;
    addRequired(p, 'LandRegion');
    addRequired(p, 'Lonlim');
    addRequired(p, 'Latlim');
    addOptional(p, 'LongitudeOrigin', 0, @isnumeric);
    parse(p, landRegion, lonlim, latlim, varargin{:});

    landRegion = p.Results.LandRegion;
    lonlim = p.Results.Lonlim;
    latlim = p.Results.Latlim;
    lonOrigin = p.Results.LongitudeOrigin;

    if ~isa(landRegion, 'polyshape')
        isPolygon = false;
        landRegion = polyshape(landRegion);
    else
        isPolygon = true;
    end

    subtractingRegionLat = [latlim(2); latlim(2); latlim(1); latlim(1)];
    subtractingRegionLon = [lonlim(1); lonlim(2); lonlim(2); lonlim(1)];
    [subtractingRegionLon, subtractingRegionLat] = ...
        poly2cw(subtractingRegionLon, subtractingRegionLat);
    [subtractingRegionLat, subtractingRegionLon] = ...
        flatearthpoly( ...
        subtractingRegionLat, subtractingRegionLon, ...
        lonOrigin);

    subtractingRegion = polyshape( ...
        subtractingRegionLon, subtractingRegionLat);

    landRegion = subtract(landRegion, subtractingRegion);

    if ~isPolygon
        landRegion = landRegion.Vertices;
        landRegion = closecoastline(landRegion);
    end

end
