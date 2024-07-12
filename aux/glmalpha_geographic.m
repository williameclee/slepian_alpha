%% GLMALPHA_GEOGRAPHIC
% Is an auxiliary function separated from GLMALPHA. 
% It finds G, V, and N for the given region.
%
% Authored by:
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12

function [G, V, N] = ...
        glmalpha_geographic(maxL, TH, sord, anti, rotb, ...
        ldim, bp, lp, EL, EM, xver, ...
        beQuiet, mesg)
    % Calculates the localization kernel for this domain
    % See if we can run this calculation in parallel
    canRunParallel = license('test', 'distrib_computing_toolbox') && (verLessThan('matlab', '8.2') && parpool('processes').NumWorkers > 0 ...
        || ~verLessThan('matlab', '8.2'));

    if canRunParallel

        if ~beQuiet
            disp('Running KERNELCP (parallel)');
        end

        Klmlmp = kernelcp_new(maxL, TH, sord);
    else

        if ~beQuiet
            disp('No open matlabpool. Running KERNELC (non-parallel).');
        end

        Klmlmp = kernelc(maxL, TH, sord);
    end

    if anti == 1
        % Get the complimentary region
        Klmlmp = eye(size(Klmlmp)) - Klmlmp;
    end

    if bp
        % So far we only have wanted to remove small portions of the
        % kernel.  Therefore at the moment here, we make the whole thing,
        % and the apply the bp after the fact.  However, in the future if
        % you want a bp at large L, we should modify kernelcp to only make
        % a partial kernel.

        % Remove the beginning section of the kernel
        rem = bp * L(1) ^ 2;
        Klmlmp(:, 1:rem) = [];
        Klmlmp(1:rem, :) = [];
    end

    % Calculates the eigenfunctions/values for this localization problem
    if lp
        [G, V] = eig(Klmlmp);
    elseif bp
        [Gbp, V] = eig(Klmlmp);
        G(L(1) ^ 2 + 1:end, :) = Gbp;
    end

    [V, isrt] = sort(sum(real(V), 1));
    V = fliplr(V);
    V = V(:);
    G = G(:, fliplr(isrt));

    [~, ~, ~, ~, ~, ~, ems, els, R1, R2] = addmon(maxL);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G = G(R1, :);
    % Check indexing
    difer(els(R1) - EL, [], [], mesg)
    difer(ems(R1) - EM, [], [], mesg)

    % Calculate Shannon number and compare this with the theory
    N = sum(V);

    if lp
        % Is the Shannon number right? Need the area of the region
        difer(ldim * Klmlmp(1) - N, [], [], mesg)
    elseif bp
        difer(ldim * spharea(TH) - N, [], [], mesg)
    end

    % Check if the expansion of a basis function is indeed either 1 or 0
    if lp && xver == 1
        disp('Excessive verification')
        % Is the area right? Don't be too demanding
        difer(Klmlmp(1) - abs(anti - spharea(TH)), 4, [], mesg)
        % This is a bit double up... but it's only for excessive verification
        [V1, C, ~, ~, ~, ~, GG] = localization(L, TH, sord);
        difer(V(:) - V1(:), [], [], mesg)
        % A test by expansion and orthogonality
        for index = 1:length(C)
            salpha = G' * C{index}(R2);
            % Only one of these functions should get "hit"
            difer(sum(abs(salpha) > 1e-8) - 1, [], [], mesg)
        end

        % Yet another way, see LOCALIZATION
        [~, ~, ~, ~, ~, ~, ~, ~, R1, ~] = addmon(L);
        % It's only a matter of indexing and ordering
        difer(GG(R1, :) - G)
    end

    % Lets check if we need to do a rotation. The function for your
    % coordinates should have this functionality if it's needed.
    try
        rotb = feval(dom, 'rotated');
    catch
    end

    % Now do the rotation
    if length(rotb) == 1 && rotb
        % Get the rotation parameters to rotate G. Note, the region
        % rotation angles that we return from the functions (lonc, latc)
        % are the same regardless of if we did a buffer, as they pertain
        % to the original region
        [~, lonc, latc] = eval(sprintf('%s()', dom));
        G = rotateGp(G, lonc, latc);
    end

end
