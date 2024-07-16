% [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=GLMALPHA(TH,L,sord,blox,upco,resc,J,anti)
%
% Returns an (lm)X(alpha) matrix with unit-normalized spherical harmonic
% coefficients of the BANDLIMITED or PASSBAND Slepian functions of the
% SINGLE or DOUBLE polar cap, or of a GEOGRAPHICAL region of
% interest. Only in the geographical case are the eigenvalues automatically
% sorted; if not, the column dimension is always block-ordered by virtue of the
% construction. The matrix G is orthogonal, G'*G is the identity. In column
% truncated form, G(:,1:J)*G(:,1:J)' is not the identity also, but rather
% a projection with eigenvalues 0 and 1.
%
% Should put an option to save only the essentials up to a certain truncation
%
% INPUT:
%
% TH       Angular extent of the spherical cap, in degrees OR
%          'england', 'eurasia',  'namerica', 'australia', 'greenland',
%          'africa', 'samerica', 'amazon', 'orinoco', 'antarctica',
%          'contshelves', 'alloceans',
%          OR: [lon lat] an ordered list defining a closed curve [degrees],
%          OR: {'region' buf} where buf is the distance in degrees that
%          the region outline will be enlarged by BUFFERM
% L        Bandwidth (maximum angular degree), or passband (two degrees)
% sord     1 Single polar cap of radius TH [default]
%          2 Double polar cap, each of radius TH
%          N Splining smoothness for geographical regions [default: 10]
%
% The following options are only for axisymmetric polar caps:
%
% blox     0 Standard (lm) row ordering, l=0:L, m=-l:l as ADDMOUT [default]
%          1 Block-diagonal row ordering, m=[0 -1 1 -2 2 ... -L L], l=abs(m):L
% upco     +ve fraction of unit radius for upward continuation [default: 0]
%          -ve fraction of unit radius for downward continuation
% resc     0 Not rescaled [default]
%          1 Rescaled to have a unit integral over the unit sphere
%            (this is only relevant for the down/upward continued functions)
%
% Though there are three more options that hold only for the geographic cases
%
% J        The number of eigenfunctions that are being asked (and saved);
% anti     1 get the opposite of the region you specify
%          0 get exactly the region that you specify [default]
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients; note how
%          LOCALIZATION delivers these as LMCOSI arrays into PLM2XYZ
% V        The eigenvalues in this ordering (not automatically sorted)
% EL       The vector with spherical harmonic degrees as first index of G
% EM       The vector with spherical harmonic orders as first index of G
% N        The Shannon number
% GM2AL    The sum over all orders of the squared coefficients, i.e. the
%          TOTAL power, NOT the power spectral density
% MTAP     The order of the eigentapers, if the region is axisymmetric
% IMTAP    The rank within that particular order of the eigentapers
%
% EXAMPLE:
%
% glmalpha('demo1') % Illustrates the block sorting and checks unitarity
% glmalpha('demo2') % Makes a coupling kernel a la BCOUPLING
% glmalpha('demo3') % Calculates something and uses PLOTSLEP to plot
%
% SEE ALSO:
%
% GLMALPHAPTO, ADDMOUT, ADDMON, KERNELC, GALPHA, DLMLMP, GLM2LMCOSI, LOCALIZATION
%
% Region functions such as ANTARCTICA have a default behavior to indicate if
% their eigenfunctions should be rotated (e.g. back to a pole). If you want
% eigenfunctions for the region at the equator then rotate them back after
% the fact using ROTATEGP.
%
% Last modified by plattner-at-alumni.ethz.ch, 10/11/2016
% Last modified charig-at-princeton.edu, 06/27/2016
% Last modified by fjsimons-at-alum.mit.edu, 12/01/2017

% Should be able to update this to retain the rank order per m as well as
% the global ordering. Does this work for the whole-sphere? In that case,
% should really want G to be the identity - all though of course,
% anything works, too. You don't get necessarily the spherical harmonics
% back...
function varargout = glmalpha_new(varargin)
    % Add path to the auxiliary functions
    addpath(fullfile(fileparts(mfilename('fullpath')), 'aux'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'demos'));

    % Parse inputs
    [TH, L, sord, blox, upco, resc, J, rotb, anti, dom, ...
         forceNew, saveData, beQuiet] = ...
        parseinputs(varargin);

    %% Demos
    if strcmp(TH, 'demo1')
        glmalpha_demo1;
        return
    elseif strcmp(TH, 'demo2')
        glmalpha_demo2;
        return
    elseif strcmp(TH, 'demo3')
        glmalpha_demo3;
        return
    end

    %% Initialisation
    defval('mesg', 'GLMALPHA Check passed')
    % Hold all messages
    mesg = NaN;

    % Figure out if it's lowpass or bandpass
    lp = isscalar(L);
    bp = length(L) == 2;

    if ~(lp || bp)
        error('The degree range is either one or two numbers')
    end

    maxL = max(L);

    % The spherical harmonic dimension
    ldim = (L(2 - lp) + 1) ^ 2 - bp * L(1) ^ 2;
    defval('J', ldim)

    % Output file
    [outputPath, outputFileExists, GM2AL, MTAP, IMTAP, xver, ~] = ...
        getoutputfile(TH, L, sord, blox, upco, resc, J, anti, ...
        lp, bp, ldim, [], [], [], [], dom);

    if outputFileExists && ~forceNew
        warning('off', 'MATLAB:load:variableNotFound')
        load(outputPath, 'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')

        defval('GM2AL', []);
        defval('MTAP', []);
        defval('IMTAP', []);

        if ~beQuiet
            fprintf('%s loading %s\n', upper(mfilename), outputPath)
        end

        varargout = {G, V, EL, EM, N, GM2AL, MTAP, IMTAP};

        return

    end

    %% Main
    % Find row indices into G belonging to the orders
    [EM, EL, ~, blkm] = addmout(maxL);

    % % Find increasing column index; that's how many belong to this order
    % % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    % % The middle bit is twice for every nonzero order missing
    % % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
    % %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
    % %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
    % % This should be the same for L and [0 L]
    % alpha = cumsum( ...
    %     [1 L(2 - lp) - bp * L(1) + 1 ...
    %      gamini(L(2 - lp) - bp * (L(1) - 1), bp * 2 * L(1)) ...
    %      gamini(L(2 - lp) - bp * L(1):-1:1, 2)]);

    % For GEOGRAPHICAL REGIONS or XY REGIONS
    if isstring(TH) || ischar(TH) || length(TH) > 1
        % Find the function in the aux/ subdirectory
        [G, V, N] = ...
            glmalpha_geographic(maxL, TH, sord, anti, rotb, ...
            ldim, bp, lp, EL, EM, xver, ...
            beQuiet, forceNew, mesg);

        G = G(:, 1:J);
        V = V(1:J);

        try
            save(outputPath, '-v7.3', 'G', 'V', 'EL', 'EM', 'N')
        catch
            save(outputPath, 'G', 'V', 'EL', 'EM', 'N')
        end

    else
        % For AXISYMMETRIC REGIONS
        [G, V, EL, EM, N, GM2AL, MTAP, IMTAP] = ...
            glmalpha_axisymmetric(TH, sord, L, lp, bp, EM, EL, blkm, ...
            blox, upco, xver);

        if ~strcmp(outputPath, 'neveravailable') && saveData
            % Save the results if it isn't a geographical region
            % If the variable is HUGE you must use the -v7.3 flag, if not, you
            % can safely omit it and get more backwards compatibility
            try
                save(outputPath, '-v7.3', 'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')
            catch
                save(outputPath, 'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')
            end

        end

    end

    %% Returning requested variables
    varargout = {G, V, EL, EM, N, GM2AL, MTAP, IMTAP};

end

%% Subfunctions
function varargout = parseinputs(vArargin)
    THDefault = 30;
    LDefault = 18;
    sordDefault = [];
    bloxDefault = 0;
    upcoDefault = 0;
    rescDefault = 0;
    JDefault = [];
    rotbDefault = 0;
    antiDefault = 0;

    p = inputParser;
    addOptional(p, 'TH', THDefault, ...
        @(x) ischar(x) || iscell(x) || isscalar(x) || isnumeric(x) || isempty(x));
    addOptional(p, 'L', LDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'sord', sordDefault, ...
        @(x) isnumeric(x) || iscell(x) || isempty(x));
    addOptional(p, 'blox', bloxDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'upco', upcoDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'resc', rescDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'J', JDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'rotb', rotbDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'anti', antiDefault, ...
        @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'ForceNew', false, @islogical);
    addParameter(p, 'Save', true, @islogical);
    addParameter(p, 'BeQuiet', false, @islogical);
    parse(p, vArargin{:});

    TH = conddefval(p.Results.TH, THDefault);
    L = conddefval(p.Results.L, LDefault);
    sord = conddefval(p.Results.sord, sordDefault);
    blox = conddefval(p.Results.blox, bloxDefault);
    upco = conddefval(p.Results.upco, upcoDefault);
    resc = conddefval(p.Results.resc, rescDefault);
    J = conddefval(p.Results.J, JDefault);
    rotb = conddefval(p.Results.rotb, rotbDefault);
    anti = conddefval(p.Results.anti, antiDefault);
    forceNew = p.Results.ForceNew;
    saveData = p.Results.Save;
    beQuiet = p.Results.BeQuiet;

    if any(strcmp(TH, {'anti'}))
        anti = 1;
    end

    dom = [];

    varargout = ...
        {TH, L, sord, blox, upco, resc, J, rotb, anti, dom, ...
         forceNew, saveData, beQuiet};
end

function [outputPath, outputFileExists, GM2AL, MTAP, IMTAP, xver, dom] = ...
        getoutputfile(TH, L, sord, blox, upco, resc, J, anti, lp, bp, ldim, GM2AL, MTAP, IMTAP, xver, dom)
    % Just get the file names here
    if upco == 0 && resc == 0

        if ~(isstring(TH) || ischar(TH) || iscell(TH)) && isscalar(TH)
            domainType = 'cap';
        elseif isstring(TH) || ischar(TH) || iscell(TH)
            domainType = 'geographic';
        else
            domainType = 'lonlat';
        end

        outputFolder = fullfile(getenv('IFILES'), 'GLMALPHA');
        defval('xver', 0)

    else
        domainType = 'pto';

        outputFolder = fullfile(getenv('IFILES'), 'GLMALPHAPTO');
        defval('xver', 1)
    end

    switch domainType
        case 'cap'
            defval('sord', 1)

            if lp
                outputFile = sprintf('glmalpha-%i-%i-%i-%i.mat', ...
                    TH, L, sord, blox);
            elseif bp
                outputFile = sprintf('glmalphabl-%i-%i-%i-%i-%i.mat', ...
                    TH, L(1), L(2), sord, blox);
            end

            % Initialize ordering matrices
            MTAP = zeros([1, ldim]);
            IMTAP = zeros([1, ldim]);
        case 'geographic'
            defval('sord', 10)
            defval('buf', 0)
            defval('J', ldim)

            if iscell(TH) && length(TH) >= 3
                regionString = [TH{1}, dataattrchar("Buffer", TH{3:end})];
                dom = TH{1};
            else

                if iscell(TH) && length(TH) == 2
                    dom = TH{1};
                    buf = TH{2};
                else
                    dom = TH;
                    buf = 0;
                end

                if iscell(sord)
                    sord = sord(~cellfun('isempty', sord));
                    regionString = [dom, dataattrchar('Buffer', buf, sord{:})];
                else
                    regionString = [dom, dataattrchar('Buffer', buf, 'Upscale', sord)];
                end

            end

            % end

            if lp
                outputFile = sprintf('glmalpha-%s-%i-%i.mat', ...
                    regionString, L, J);
            elseif bp
                outputFile = sprintf('glmalphabl-%s-%i-%i-%i.mat', ...
                    regionString, L(1), L(2), J);
            end

        case 'lonlat'
            defval('sord', 10)
            defval('buf', 0)
            defval('J', ldim)
            regionString = hash(TH, 'sha1');

            if lp
                outputFile = sprintf('glmalpha-%s-%i-%i.mat', ...
                    regionString, L, J);
            elseif bp
                outputFile = sprintf('glmalphabl-%s-%i-%i-%i.mat', ...
                    regionString, L(1), L(2), J);
            end

        case 'pto'
            % Make a hash, who cares if it's human-readable?
            outputFile = sprintf('%s.mat', hash([TH L phi theta omega], 'sha1'));
            % For excessive verification of the upco'd case
    end

    outputPath = fullfile(outputFolder, outputFile);

    if anti == 1
        % Update the file name to reflect the complimentarity of the region
        outputPath = sprintf('%s-1.mat', pref(outputPath));
    end

    outputFileExists = exist(outputPath, 'file') == 2;

end
