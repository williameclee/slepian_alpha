%% AFRICA
% Returns the longitude-latitude coordinates of the African continent,
% potentially buffered by some amount.
%
% Syntax
%   lonlat = africa(upscale, buf)
%   africa(__)
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
%   williameclee-at-arizona.edu, 07/30/2024
%   fjsimons-at-alum.mit.edu, 11/23/2011
%   charig-at-princeton.edu, 11/23/2011

function lonlat = africa(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'africa';
    c11 = [0, 37.5; 342.5, 36.5];
    cmn = [52, -35; 360, 3.5];
    xunt = [261:307, 110:259, 339:399];
    ofs = [360, 0];

    % Find/load/save the coordinates
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf, ofs);

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of Africa (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
