%% ENGLAND
% Returns the longitude-latitude coordinates of the UK, potentially 
% buffered by some amount.
% Note that Ireland and North Ireland are not included.
%
% Syntax
%   lonlat = england(upscale, buf)
%   england(__)
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
% Note
%   This region's shapce is self-intersecting, may resolve this in the 
%   future
%
% Last modified by
%   williameclee-at-arizona.edu, 07/30/2024
%   fjsimons-at-alum.mit.edu, 11/23/2011
%   charig-at-princeton.edu, 11/23/2011

function lonlat = england(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'england';
    c11 = [353, 59; 0, 59];
    cmn = [360, 49.75; 2, 49.75];
    xunt = [61:68, 9:60, 72:80];
    ofs = [0, 360];

    % Find/load/save the coordinates
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf, ofs);

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of England (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
