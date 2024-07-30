%% GREENLAND
% Returns the longitude-latitude coordinates of greenland,
% potentially buffered by some amount.
% Note that the nearby island of Ellesmere is not included and may be
% subtracted.
%
% Syntax
%   lonlat = greenland(upscale, buf)
%   lonlat = greenland(upscale, buf, nearby)
%   greenland(__)
%
% Inputs
%   upscale - The times of spline-upscaling applied to the coordinates
%       The default value is 0 (no upscaling)
%   buf - The size of the buffer from the coastline in degrees
%       The value can be positive (buffering outwards) or negative
%       (buffering inwards)
%       The default value is 0 (no buffer)
%  nearby - Whether to subtract the nearby island of Ellesmere
%       The default value is false
%
% Outputs
%   lonlat - Closed-curved coordinates of the continent
%       The coordinates are in the form of [longitude(:), latitude(:)] in
%       degrees
%
% Note
%   The nearby functionality is not properly tested, but it should work
%   fine (williameclee-at-arizona.edu, 07/30/2024)
%
% Last modified by
%   williameclee-at-arizona.edu, 07/30/2024
%   charig-at-princeton.edu, 07/10/2014
%   fjsimons-at-alum.mit.edu, 11/23/2011

function lonlat = greenland(varargin)
    % Parse inputs
    [upscale, buf, nearby] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'greenland';
    c11 = [286.7, 83.75];
    cmn = [349, 59.75];
    xunt = 1:352;

    % Do it! Make it, load it, save it
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf);

    % Subtract the nearby island of Ellsemere
    if nearby
        % NOTE: if you don't have the nearby stuff made at the desired 
        % buffer then this will just end in a recursive error. i.e. You 
        % should make Ellesmere first before trying to subtract it.
        lonlat2 = ellesmereg(10, buf);
        [lon, lat] = polybool('subtraction', ...
            lonlat(:, 1), lonlat(:, 2), lonlat2(:, 1), lonlat2(:, 2));
        lonlat = [lon, lat];

        if buf >= 2
            lonlat2 = baffing(10, 1.5);
            [lon, lat] = polybool('subtraction', ...
                lonlat(:, 1), lonlat(:, 2), lonlat2(:, 1), lonlat2(:, 2));
            lonlat = [lon, lat];
        end

    end

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of Greenland (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
