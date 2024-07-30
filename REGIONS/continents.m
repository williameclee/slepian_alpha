%% CONTINENTS
% Returns the longitude-latitude coordinates of all major continents,
% potentially buffered by some amount.
% Note that the Antarctica is not included.
%
% Syntax
%   lonlat = continents(upscale, buf)
%   continents(__)
%
% Inputs
%   upscale - The times of spline-upscaling applied to the coordinates
%       The default value is 0 (no upscaling)
%   buf - The size of the buffer from the coastline in degrees
%       The value can be positive (buffering outwards) or negative
%       (buffering inwards)
%       The default value is 0 (no buffer)
%
% Outputs
%   lonlat - Closed-curved coordinates of the continent
%       The coordinates are in the form of [longitude(:), latitude(:)] in
%       degrees
%
% Last modified by
%   williameclee-at-arizona.edux, 07/30/2024

function lonlat = continents(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);

    %% Combining the regions
    continentList = {'namerica', 'samerica', 'africa', 'eurasia', 'australia', 'greenland'};

    for iCont = 1:length(continentList)
        lonlatc = feval(continentList{iCont}, upscale, buf);

        if strcmp(continentList{iCont}, 'australia')
            lonlatc = [lonlatc(:, 1) + 360, lonlatc(:, 2)];
        end

        if iCont == 1
            lonlat = lonlatc;
        else
            lonlat = [lonlat; NaN, NaN; lonlatc]; %#ok<AGROW>
        end

    end

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of major landmasses (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
