%% AUSTRALIA
% Returns the longitude-latitude coordinates of Australia, potentially
% buffered by some amount.
% Note that the Tasmania is not included.
%
% Syntax
%   lonlat = australia(upscale, buf)
%   australia(__)
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

function lonlat = australia(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'australia';
    c11 = [112, -10.5];
    cmn = [154, -39];
    xunt = 2:166;

    % Find/load/save the coordinates
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf);

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of Australia (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
