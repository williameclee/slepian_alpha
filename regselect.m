% XY=REGSELECT(regn,c11,cmn,xunt,res,buf,ofs)
%
% Returns coordinates of a certain named specified region, and will
% help you make an appropriate (buffered) selection before saving it.
%
% INPUT:
%
% regn     The region name string, e.g. 'greenland'
% c11      The (lon,lat) coordinates of the top left corner
%          or a [lon(:)] set of coordinates for a region
% cmn      The (lon,lat) coordinates of the bottom right corner
%          or a [lat(:)] set of coordinates for a region
% xunt     The particuluar indices required for this region
% upscale  0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
% ofs      Offset(s) to be considered by PLOTCONT [default: none]
%
% OUTPUT:
%
% XY       The requested coordinates
%
% Last modified by charig-at-email.arizona.edu, 12/05/2017
% Last modified by maxvonhippel@email.arizona.edu, 12/05/2017
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2022

function XY = regselect(regn, c11, cmn, varargin)
    %% Parsing inputs
    p = inputParser;
    p.addRequired('regn', @ischar)
    p.addRequired('c11', @isnumeric)
    p.addRequired('cmn', @isnumeric)
    p.addOptional('xunt', [], @isnumeric)
    p.addOptional('upscale', 0, @(x) isnumeric(x) || isempty(x))
    p.addOptional('buf', 0, @(x) isnumeric(x) || isempty(x))
    p.addOptional('ofs', 0, @(x) isnumeric(x) || isempty(x))
    p.addParameter('spd', 0.01, @isnumeric)
    parse(p, regn, c11, cmn, varargin{:})

    Regn = capitalise(p.Results.regn);
    c11 = p.Results.c11;
    cmn = p.Results.cmn;
    xunt = p.Results.xunt;
    upscale = p.Results.upscale;

    if isempty(upscale)
        upscale = 0;
    end

    buf = p.Results.buf;

    if isempty(buf)
        buf = 0;
    end

    ofs = p.Results.ofs;
    if isempty(ofs)
        ofs = 0;
    end
    spd = p.Results.spd;

    % The directory where you keep the coordinates
    dataFolder = fullfile(getenv('IFILES'), 'COASTS');

    %% Look for the data file
    % Revert to original name if unbuffered
    if upscale == 0 && buf == 0
        dataFile = sprintf('%s.mat', Regn);
    elseif buf == 0
        dataFile = sprintf('%s-%i.mat', Regn, upscale);
    else
        dataFile = sprintf('%s-%i-%g.mat', Regn, upscale, buf);
    end

    dataFile = fullfile(dataFolder, dataFile);

    % If you already have a file
    if exist(dataFile, 'file') == 2
        load(dataFile, 'XY')
        % fprintf('%s loading %s\n', upper(mfilename), dataFile)
        return
    end

    %%
    if upscale ~= 0
        XY = bezier(eval(sprintf('%s(0)', regn)), upscale);
    else
        figure
        clf
        % This could be 1x2, 2x1 or 2x2 or 2x3 indeed
        if length(c11) <= 3 && length(cmn) <= 3
            XY = [];

            for index = 1:size(c11, 1)
                [~, handl, XY2] = plotcont(c11(index, :), cmn(index, :), [], ofs(index));
                % Get rid of common NaNs
                XY2 = XY2(~isnan(XY2(:, 1)) & ~isnan(XY2(:, 2)), :);
                delete(handl)
                XY = [XY; XY2]; %#ok<AGROW>
                clear XY2
            end

        else
            XY(:, 1) = c11;
            XY(:, 2) = cmn;
        end

        % Get rid of common NaNs
        XY = XY(~isnan(XY(:, 1)) & ~isnan(XY(:, 2)), :);

        % If this is antarctica, we want to rotate to the equator
        if strcmp(regn, 'antarctica')
            XY = rotateantarctica(XY);
        end

        % Now make sure the distances aren't huge
        [X, Y] = penlift(XY(:, 1), XY(:, 2), 3);
        XY = [X, Y];

        % Some handiwork specifically for antarctica
        if strcmp(regn, 'antarctica')
            XY = handiworkantarctica(XY);
        end

        % Experiment with cutoff -> See "check this out"
        if ~isempty(xunt)
            XY = XY(xunt, :);
        end

        % Definitely close the contour again
        XY = closecoastline(XY);

        % For Ellesmere Island we want a three part polygon
        if strcmp(regn, 'ellesmere')
            XY = handiworkellesmere(XY);
        end

        % And definitely make this go clockwise
        [X2, Y2] = poly2cw(XY(:, 1), XY(:, 2));
        XY = [X2, Y2];
        clear X2 Y2

        plot(XY(:, 1), XY(:, 2), 'LineW', 2, 'Color', 'k');
        axis equal

        hold on
        curvecheck(XY(:, 1), XY(:, 2), spd);
        hold off
    end

    %% Buffer XY
    if buf ~= 0

        if buf > 0
            % Here we assume you now have a 2015 or later version of MATLAB
            % so we can interpret positive buffers as outPlusInterior. In
            % prior versions, this used to be just 'out' and require a fix.
            inout = 'outPlusInterior';
        else
            inout = 'in';
        end

        disp('Buffering the coastlines... this may take a while');
        [latBuf, lonBuf] = bufferm(XY(:, 2), XY(:, 1), abs(buf), inout);
        % You might look into REDUCEM to make this easier
        % XY=[LonB{2}+180*[LonB{2}<0] LatB{2}];

        % Note that, if due to BEZIER there might be a pinched-off loop in
        % the XY you will get an extra NaN and will need to watch it
        % If this shouldn't work, keep it all unchanged in other words

        % Periodize the right way
        XY = [lonBuf + 360 * any(lonBuf < 0), latBuf];
    end

    %% Save data
    save(dataFile, 'XY')

end

%% Subfunctions
function XYRotated = rotateantarctica(XY)
    % Find the geographical center and the area
    [~, centerLat] = rcenter([XY(:, 1), XY(:, 2)]);
    % Just make it not straddle the equator; flip the sign;
    % stay close to the value calculated by RCENTER
    centerLon = -45;
    % Convert to Cartesian coordinates
    [X, Y, Z] = sph2cart(deg2rad(XY(:, 1)), deg2rad(XY(:, 2)), 1);
    % Apply the rotation to put it on or near the equator
    XYZRotated = (rotz(deg2rad(centerLon)) * ...
        roty(deg2rad(-centerLat)) * ...
        [X(:), Y(:), Z(:)]')';
    % See LOCALIZATION and KLMLMP2ROT for the counterrotation
    % Transform back to spherical coordinates
    [lonRotated, latRotated] = cart2sph( ...
        XYZRotated(:, 1), XYZRotated(:, 2), XYZRotated(:, 3));
    XYRotated = [rad2deg(lonRotated), rad2deg(latRotated)];
    % Here, lonc and latc are the same values given out from antarctica
    % for the rotation back to the pole.
end

function XY = handiworkantarctica(XY)
    XY = [flipud(XY(p + 2:p + 3, :));
          XY([[(1:p(1))]';
         (p(1) + 5:size(XY, 1))'], :)];
    xx = XY(:, 1);
    yy = XY(:, 2);
    d = sqrt((xx(2:end) - xx(1:end - 1)) .^ 2 + (yy(2:end) - yy(1:end - 1)) .^ 2);
    XY = XY(d > 1e-14, :);
    XY = [XY(2:227, :); XY(232:end, :); XY(2, :)];
end

function XY = handiworkellesmere(XY)
    % Undo what we just did matching the end to the start
    % and insert NaNs in order to make this a 3 parts
    XY = XY(1:end - 1, :);
    XY = insert(XY, [NaN, NaN, NaN, NaN], [24, 150, 209, 335]);
    XY = reshape(XY, [], 2);
end
