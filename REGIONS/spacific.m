%% SPACIFIC
% XY = SPACIFIC(upscale, inclang)
% SPACIFIC(...) % Only makes a plot
%
% Finds the coordinates of some of the worlds' oceans so we can combine
% its localization kernel with those for the missing continents to turn
% into Slepian eigenfunctions for all of the world's oceans.
%
% INPUT:
%
% upscale  0 The standard, default values
%          N Splined values at N times the upscaleolution
% inclang  Inclination angle, determines size of polar caps
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the Atlantic ocean
%
% Last modified by williameclee-at-arizona.edu, Jul 2nd, 2024

function varargout = spacific(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, inclang, buf, moreBuf, forceReload, lonOrigin] = ...
        parsecoastinputs(varargin, 'DefaultLongitudeOrigin', 180);
    oceanTitle = 'North Pacific Ocean';
    oceanParts = ...
        {'South Pacific Ocean, eastern part', ...
     'South Pacific Ocean, western part'};

    %% Check if the data file exists
    [dataFile, ~, dataFileExists] = coastfilename(mfilename, ...
        'Upscale', upscale, 'Inclang', inclang, ...
        'Buffer', buf, 'MoreBuffer', moreBuf);

    if dataFileExists && ~forceReload
        load(dataFile, 'XY')

        varargout = returncoastoutput(nargout, XY, oceanTitle);

        return
    end

    %% Compute the data
    [oceanPoly, oceanLatlim, oceanLonlim] = ...
        findoceanboundary(oceanParts, inclang, lonOrigin);

    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Convert the data to the right format
    [X, Y] = poly2cw(XY(:, 1), XY(:, 2));
    % [Y, X] = flatearthpoly(Y, X, lonOrigin);
    XY = removeduplicatevertices([X, Y]);

    %% Save and return required data
    save(dataFile, 'XY', '-v7.3')
    fprintf('%s saving %s\n', upper(mfilename), dataFile)

    varargout = returncoastoutput(nargout, XY, oceanTitle);

end

function coastPoly = manualadjustment(coastPoly, buf, lonOrigin)

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-40, -39], 'Lonlim', [143.5, 145.5], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 2
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-1, 0], 'Lonlim', [135, 137], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 4
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-2, 0], 'Lonlim', [145, 152], ...
            'LongitudeOrigin', lonOrigin);
    end

end
