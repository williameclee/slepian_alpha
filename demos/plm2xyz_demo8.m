%% PLM2XYZ_DEMO8
% This is demo 8 from PLM2XYZ.
%
% See also
%   PLM2XYZ
% Authored by
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-25

function output = plm2xyz_demo8(meshSize)
    v = fralmanac('EGM2008_ZeroTide', 'SHM');
    % Note that gravity does not start at zero
    % Geoid = 3 % Free-air gravity anomaly = 2
    v = plm2pot(v(1:addmup(720) - addmup(v(1) - 1), :), [], [], [], 3);
    [r, lon, lat, Plm] = plm2xyz(v);
    plotplm(r)
    % Prepare output
    output = {r, lon, lat, Plm, meshSize};
end
