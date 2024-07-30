% [Bl,dels,r,lon,lat,lmcosi]=GEOBOXCAP(L,dom,N,meshSize,act,lonc,latc)
%
% Returns the spherical-harmonic POWER SPECTRUM in the UNIT-NORMALIZED
% spherical-harmonic basis of an all-or-nothing BOXCAR over a particular
% geographical region. It does that by calculating the power spectrum in
% the 4-pi-normalized basis and then adjusting the coefficients to pick up
% the power so that it matches that in the unit-normalized basis.
%
% INPUT:
%
% L          Bandwidth of the output spectrum
% dom        'africa',  'eurasia',  'namerica', 'australia', 'greenland',
%            'samerica', 'amazon', 'orinoco', 'gpsnamerica', 'antarctica'
%            OR: [lon lat] an ordered list defining a closed curve [degrees]
% N          The splining smoothness of the geographical boundary [default: 10]
% meshSize     The resolution of the underlying spatial grid [default: Nyquist]
% act        1 Actually perform the calculations [default]
%            0 Don't, but simply return the (rotated) mask function
% lonc,latc  Rotate coordinates by this amount
%
% OUTPUT:
%
% Bl         The power spectrum
% dels       The spherical harmonic degrees
% r          The "mask" with the "continent function"
% lon,lat    The grid on which  this is defined
% lmcosi     The complete set of spherical harmonic expansion coefficients
%
% EXAMPLE:
%
% [Bl1,dels1]=geoboxcap(18,'australia',[],[]);
% [Bl2,dels2]=geoboxcap(18,'australia',[],1);
%
% SEE ALSO:
%
% BPBOXCAP, KERNELC, GAMMAP
%
% Last modified by fjsimons-at-alum.mit.edu, 05/13/2013
function varargout = geoboxcap(varargin)

    %% Parsing inputs and initialisation
    c11cmn = [0, 90, 360, -90]; % default to global

    % Parse inputs
    p = inputParser;
    addOptional(p, 'L', 18);
    addOptional(p, 'Domain', 'australia', ...
        @(x) (ischar(x) || iscell(x)));
    addOptional(p, 'Upscale', 1);
    addOptional(p, 'MeshSize', -1);
    addOptional(p, 'act', 1);
    addOptional(p, 'lonc', 0);
    addOptional(p, 'latc', 0);
    addOptional(p, 'Inverse', false, ...
        @islogical);
    addOptional(p, 'ForceNew', false, ...
        @islogical);
    parse(p, varargin{:});

    if ~isempty(p.Results.L)
        L = p.Results.L;
    else
        L = 18;
    end

    if ~isempty(p.Results.Domain)
        dom = p.Results.Domain;
    else
        dom = 'australia';
    end

    if ~isempty(p.Results.Upscale)
        upscale = p.Results.Upscale;
    else
        upscale = 1;
    end

    if ~isempty(p.Results.MeshSize)
        meshSize = p.Results.MeshSize;
    else
        meshSize = -1;
    end

    if ~isempty(p.Results.act)
        act = p.Results.act;
    else
        act = 1;
    end

    if ~isempty(p.Results.lonc)
        lonc = p.Results.lonc;
    else
        lonc = 0;
    end

    if ~isempty(p.Results.latc)
        latc = p.Results.latc;
    else
        latc = 0;
    end

    if ~isempty(p.Results.Inverse)
        isInverted = p.Results.Inverse;
    else
        isInverted = false;
    end

    if ~isempty(p.Results.ForceNew)
        forceNew = p.Results.ForceNew;
    else
        forceNew = false;
    end

    % Make a grid of ones and zeroes depending on the desired resolution
    degN = 180 / sqrt(L * (L + 1));

    if meshSize < 0
        meshSize = degN;
    end

    % Find the coordinates of the region
    if ischar(dom)
        % Run the named function to return the coordinates
        XY = feval(dom, upscale);
    elseif iscell(dom)

        if isocean(dom{1})
            XY = feval(dom{:});
        else
            XY = feval(dom{1}, dom{3}, dom{5});
        end

    else
        % Get the coordinates as defined from the input in degrees
        XY = dom;
    end

    % Make sure the coordinates make sense
    if any(XY(:, 1) > 360) || any(XY(:, 1) < 0)
        [Y, X] = flatearthpoly(XY(:, 2), XY(:, 1), 180);
        XY = [X, Y];
    end

    %% Find data file
    [dataFile, fileExists] = getoutputfile(L, dom, upscale, meshSize, lonc, latc, isInverted);

    if fileExists && ~forceNew
        fprintf('GEOBOXCAP loading %s \n', dataFile)
        load(dataFile, 'Bl', 'dels', 'r', 'lon', 'lat', 'lmcosi')

        if nargout == 0
            plotbox(lon, lat, r, XY, lmcosi)
        else
            varargout = {Bl, dels, r, lon, lat, lmcosi};
        end

        return
    end

    %% Main code
    % The number of longitude and latitude grid points that will be computed
    nlon = ceil((c11cmn(3) - c11cmn(1)) / meshSize + 1);
    nlat = ceil((c11cmn(2) - c11cmn(4)) / meshSize + 1);
    % Lon and lat grid vectors in degrees
    lon = linspace(c11cmn(1), c11cmn(3), nlon);
    lat = linspace(c11cmn(2), c11cmn(4), nlat);
    % Make the input grid
    r = zeros(nlat, nlon);
    [LON, LAT] = meshgrid(lon, lat);

    % Now decide if we're inside or outside of the region
    % This isn't going to work for continents straddling the date line as in
    % Africa, so, rotate such cases out of the way! E.g. for the Volta basin,
    % stick on lonc=-10;
    if lonc ~= 0 || latc ~= 0
        XYor = XY;
        [X, Y, Z] = sph2cart(deg2rad(XY(:, 1)), deg2rad(XY(:, 2)), 1);
        xyzp = (rotz(deg2rad(lonc)) * roty(deg2rad(-latc)) ...
            * [X(:), Y(:), Z(:)]')';
        [phi, piminth] = cart2sph(xyzp(:, 1), xyzp(:, 2), xyzp(:, 3));
        lon = rad2deg(phi);
        lat = rad2deg(piminth);
        XY = [lon, lat];
    end

    r(inpolygon(LON, LAT, XY(:, 1), XY(:, 2))) = 1;

    % Invert the mask if requested
    if isInverted
        r = 1 - r;
    end

    if act == 1
        % And now do the spherical harmonic transform but only to L
        % This takes time!
        lmcosi = xyz2plm(r, L);

        if lonc ~= 0 || latc ~= 0
            % Now have the pleasure to rotate this back
            if latc == 0
                % Only a longitudinal rotation, much faster!
                [C, S] = rotcof(lmcosi(:, 3), lmcosi(:, 4), -lonc * pi / 180);
                lmcosi(:, 3) = C;
                lmcosi(:, 4) = S;
            else
                % A longitudinal and latitudinal rotation
                lmcosi = plm2rot(lmcosi, -lonc, latc, 0);
            end

        end

        % Calculate the power spectral density in the 4pi basis
        [Bl, dels] = plm2spec(lmcosi, 2);
    else
        [Bl, dels, lmcosi] = deal(NaN);
    end

    % Check the B0 term which should equal the area^2 divided by (4pi)^2 in the
    % 4pi-normalized basis, where the Y00 term equals 1
    A1 = 4 * pi * spharea(XY);
    A2 = 4 * pi * areaint(XY(:, 2), XY(:, 1));
    fprintf('GEOBOXCAP A: %6.3f ; SPHAREA A: %6.3f ; AREAINT A: %6.3f\n', ...
        4 * pi * sqrt(Bl(1)), A1, A2)

    % Make the adjustment so that this power spectrum is like the one from
    % BPBOXCAP, i.e. so that it is also quoted in the unit-normalized
    % harmonics.
    Bl = Bl * 4 * pi;

    % If these ever used to build loops etc they must be a row
    dels = dels(:)';

    %% Save and return requested data
    save(dataFile, 'Bl', 'dels', 'r', 'lon', 'lat', 'lmcosi')
    fprintf('%s saving %s\n', upper(mfilename), dataFile)

    if nargout == 0

        if ~exist('XYor', 'var')
            XYor = XY;
        end

        plotbox(LON, LAT, XY, XYor, lmcosi)
    else
        varns = {Bl, dels, r, lon, lat, lmcosi};
        varargout = varns(1:nargout);
    end

end

%% Subfunctions
function plotbox(LON, LAT, XY, XYor, lmcosi)
    figure
    clf
    fridplot(LON, LAT)
    xlim(xpand(minmax(LON)))
    ylim(xpand(minmax(LAT)))
    hold on
    plot(XY(:, 1), XY(:, 2), 'r')
    hold off

    figure
    clf
    rp = plm2xyz(lmcosi);
    imagefnan([0 90], [360 -90], rp)
    hold on
    % Good enough for plotting
    XYor(:, 1) = XYor(:, 1) + 360 * (XYor(:, 1) < 0);
    [X, Y] = penlift(XYor(:, 1), XYor(:, 2));
    plot(X, Y, 'r')
    hold off
end

function [dataFile, fileExists] = getoutputfile(L, dom, upscale, meshSize, lonc, latc, isInverted)
    dataFolder = fullfile(getenv('GEOBOXCAP'));

    if isempty(dataFolder)
        dataFolder = fullfile(getenv('IFILES'), 'GEOBOXCAP');
    end

    if ~exist(dataFolder, 'dir')
        warning('Data folder does not exist. Creating it now.')
        mkdir(dataFolder)
    end

    % dataFileAttrXY = cell(1, 5);

    if isInverted

        if iscell(dom)
            dataFileAttrXY = ['I', capitalise(dom{1})];
        else
            dataFileAttrXY = ['I', capitalise(dom)];
        end

    else

        if iscell(dom)
            dataFileAttrXY = capitalise(dom{1});
        else
            dataFileAttrXY = capitalise(dom);
        end

    end

    if iscell(dom)

        if length(dom) >= 3
            dataFileAttrXY = [dataFileAttrXY, dataattrchar("Buffer", dom{3:end})];
        elseif length(dom) >= 2
            dataFileAttrXY = [dataFileAttrXY, dataattrchar("Buffer", dom{2})];
        end

    end

    dataFileAttrXY = [dataFileAttrXY, '-', num2str(upscale), '-', num2str(L)];

    dataFileAttrRes = cell(1, 3);
    dataFileAttrRes{1} = num2str(meshSize);

    if lonc ~= 0 || latc ~= 0
        dataFileAttrRes{2} = sprintf('-%f', lonc);
        dataFileAttrRes{3} = sprintf('-%f', latc);
    end

    dataFile = [dataFileAttrXY, '_', dataFileAttrRes{:}, '.mat'];

    dataFile = fullfile(dataFolder, dataFile);
    fileExists = exist(dataFile, 'file') == 2;
end
