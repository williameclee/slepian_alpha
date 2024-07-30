%% NAMERICA
% Returns the longitude-latitude coordinates of the North America 
% continent, potentially buffered by some amount.
%
% Syntax
%   lonlat = namerica(upscale, buf)
%   nnamerica(__)
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

function lonlat = namerica(varargin)
    % Parse inputs
    [upscale, buf] = parsedomaininputs(varargin);
    % Parameters that make this the region in question
    domainName = 'namerica';
    c11 = [191, 74.5];
    cmn = [304.5, 12];
    % xunt = [385:717, 721:1004];
    xunt = [385:717, 722:1004];

    % Find/load/save the coordinates
    lonlat = regselect(domainName, c11, cmn, xunt, upscale, buf);

    % Plot the result if no output is requested
    if nargout > 0
        return
    end

    figure(10)
    % Specify a figure number so there won't be a new figure each time
    set(gcf, 'Name', ...
        sprintf('Coordinates of North America (%s)', upper(mfilename)))
    plot(lonlat(:, 1), lonlat(:, 2), 'k-')
    axis image
    grid on
end
