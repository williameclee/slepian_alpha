% [phint,thp,php,forreal]=PHICURVE(collon,th)
%
% Finds the longitude crossings and (thus) integration domains of a
% closed curve parameterized in colatitude/longitude space at certain
% query points of colatitude. Note that the points must be given as
% running around the curve clockwise or anticlockwise. Note that the
% region is treated as occupying flat Cartesian geometry (which is,
% admittedly, a bit of a limitation). Could be curves separated by NaNs.
%
% INPUT:
%
% collon          Colatitude/longitude of the closed curve [degrees]
% th              Colatitude at which crossings are required [degrees]
%
% OUTPUT:
%
% phint           A matrix with crossings/intervals and zeroes, of
%                 dimensions MxN where M=length(th) and N can be anything
%                 depending on the number of oscillations of the curve
% thp             Colatitude matrix for hatched plotting, if possible
% php             Longitude matrix for hatched plotting, if possible
% forreal         Indices to the ones that are for real (could be AT zero)
%
% EXAMPLE:
%
% phicurve('demo1','africa') % A geographic region
% phicurve('demo2') % A random blob
%
% SEE ALSO:
%
% DPHREGION, DPHPATCH, DPHSQUARE, CURVECHECK
%
% Last modified by charig-at-princeton.edu, 08/14/2015
% Last modified by fjsimons-at-alum.mit.edu, 11/07/2016
% Last modified by williameclee-at-arizona.edu, 11/06/2024
function varargout = phicurve(collon, th)

    %% Deal with the demos
    if strcmp(collon, 'demo1')
        rundemo1(th)
        return
    elseif strcmp(collon, 'demo2')
        rundemo2
        return
    end

    %% Main function
    % For every th, find the relevant phint
    % This is a terrible method, but by adding a very small number to the colatitude,
    % we can avoid the case where the colatitude of a coordinate is the exact same as the crossing
    % This issue seems to only arises when th is somewhat regular
    % distFromTh = collon(:, 1) + (1e-14) - th(:)';
    distFromTh = collon(:, 1) - th(:)';
    hasCrossing = diff(sign(distFromTh));
    hasCrossing(isnan(hasCrossing)) = 0;

    if sum((hasCrossing(:) == 0) - 1) == 0
        error('Specify at least one colatitude within the data range')
    elseif sum(hasCrossing(:)) ~= 0
        distFromTh = collon(:, 1) + (any(collon(:, 1) == th(:)', 2) * 1e-14) - th(:)';
        hasCrossing = diff(sign(distFromTh));
        hasCrossing(isnan(hasCrossing)) = 0;

        if sum(hasCrossing(:)) ~= 0
            error(sprintf('Cannot find pairs of crossings'))
        end

    end

    % Now it can be the one before, or after the crossing, how about
    [iXYCrossing, iColCrossing] = find(hasCrossing);
    nCrossing = length(iXYCrossing);
    % A crossing is mistakenly counted twice if the colatitude of the coordinate and the target are the same
    nDupCrossing = sum(hasCrossing(1:end - 1, :) == 1 & hasCrossing(2:end, :) == 1, 'all') + ...
        sum(hasCrossing(1:end - 1, :) == -1 & hasCrossing(2:end, :) == -1, 'all');
    nCrossing = nCrossing - nDupCrossing;

    ijMxCrossing = sub2ind(size(hasCrossing), iXYCrossing, iColCrossing);
    % This now returns the one on the negative side of the line
    iXYCrossingU = iXYCrossing + (hasCrossing(ijMxCrossing) == -2);
    iXYCrossingL = iXYCrossing + (hasCrossing(ijMxCrossing) == 2);

    % The better logic: summation should be zero
    if mod(nCrossing, 2)
        error(sprintf('Cannot find pairs of crossings'))
    end

    % Then one point was exactly hit, this is the thN or thS case
    if length(iXYCrossingU) == 2 && all(iXYCrossingU == iXYCrossingL)
        phint = collon([iXYCrossingU(2) iXYCrossingL(2)], 2);
        thp = [th th];
        php = phint;
    else

        for iCrossing = 1:nCrossing
            % In case you have a node immediate followed by a crossing
            if iXYCrossingU(iCrossing) == iXYCrossingL(iCrossing)
                phint(iCrossing) = NaN;
            else
                phint(iCrossing) = interp1(distFromTh([iXYCrossingU(iCrossing) iXYCrossingL(iCrossing)], iColCrossing(iCrossing)), ...
                    collon([iXYCrossingU(iCrossing) iXYCrossingL(iCrossing)], 2), 0, 'linear');
            end

        end

        % Debate whether this is useful or not wrt to node/crossing
        %    phint=phint(~isnan(phint));
        % ACTUALLY, IF THE NAN'S ARE NOT CONSECUTIVE PAIRS GET SPECIAL CASE

        % Now rearrange back to the number of requested points
        % But there could be points with more or less than 2 crossings
        % Maximum number of times a crossing is repeated
        [a, b] = degamini(iColCrossing);
        rowj = iColCrossing;
        iColCrossing = matranges(reshape([repmat(1, length(b), 1) b']', length(b) * 2, 1))';
        pint = repmat(NaN, length(th), max(b));
        subsi = (iColCrossing - 1) * length(th) + rowj;
        pint(subsi) = phint;

        if length(b) == length(th)
            wt = 0;
            thp = reshape(gamini(th, b), 2, length(phint) / 2);
        else
            wt = 1;
            thp = [];
        end

        % Need to sort since contour may be given in any order
        phint = sort(pint, 2);

        if wt == 0
            php = reshape(phint(subsi), 2, length(iColCrossing) / 2);
        else
            php = [];
        end

        % Make them zero so the integral doesn't do anything
        forreal = ~isnan(phint);
        phint(~forreal) = 0;
    end

    varns = {phint, thp, php, forreal};
    varargout = varns(1:nargout);

end

%% Demos
function rundemo1(th)
    defval('N', 10)
    defval('th', 'namerica')
    defval('th', 'africa')
    region = th;
    eval(sprintf('XY=%s(%i);', region, N));
    collon = [90 - XY(:, 2) XY(:, 1)]; %#ok<USENS>, XY obtained from eval
    Nth = ceil(rand * 300);
    th = linspace(min(collon(:, 1)), max(collon(:, 1)), Nth);
    [~, thp, php] = phicurve(collon, th);

    plot(php, 90 - thp, 'k-')
    hold on
    plot(collon(:, 2), 90 - collon(:, 1), 'k-')
    hold off
    axis equal; grid on
    title(sprintf('Number of crossings %i', Nth))
end

function rundemo2
    [x, y] = blob(1, 1);
    collon = [y(:) x(:)];
    Nth = ceil(rand * 300);
    th = linspace(min(collon(:, 1)), max(collon(:, 1)), Nth);
    [~, thp, php] = phicurve(collon, th);

    plot(php, 90 - thp, 'k-')
    hold on
    plot(collon(:, 2), 90 - collon(:, 1), 'k-')
    hold off
    axis equal; grid on
    title(sprintf('Number of crossings %i', Nth))
end
