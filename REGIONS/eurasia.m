%% EURASIA
% Returns the longitude-latitude coordinates of the Eurasia continent,
% potentially buffered by some amount.
% Note that the longitude is offset by 360 degrees and more to make the 
% curve continuous.
%
% Syntax
%   lonlat = eurasia(upscale, buf)
%   eurasia(__)
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

function lonlat = eurasia(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'eurasia';
    c11 = [0, 77.5; 350, 50];
    cmn = [180, 8; 360, 36];
    xunt = [420:827, 914:1023, 1556:1605, 1030:1548];
    ofs = [360, 0, 360];

    % Find/load/save the coordinates
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf, ofs);
    lonlat = closecoastline(lonlat);

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
