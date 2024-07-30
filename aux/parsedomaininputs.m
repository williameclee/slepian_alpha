function [upscale, buffer, nearby] = ...
        parsedomaininputs(inputs, varargin)
    %% Default values
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p, 'Inputs', @(x) iscell(x));
    addOptional(p, 'DefaultUpscale', 0, @(x) isnumeric(x));
    addOptional(p, 'DefaultBuffer', 0, @(x) isnumeric(x));
    % nearby is only used in greenland.m
    addOptional(p, 'DefaultNearby', false, ...
        @(x) isnumeric(x) || islogical(x));
    parse(p, inputs, varargin{:});

    inputs = p.Results.Inputs;
    defaultUpscale = p.Results.DefaultUpscale;
    defaultBuffer = p.Results.DefaultBuffer;
    defaultNearby = p.Results.DefaultNearby;

    clear p

    %% Real inputs
    p = inputParser;
    p.KeepUnmatched = true;
    addOptional(p, 'Upscale', defaultUpscale, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', defaultBuffer, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Nearby', defaultNearby, ...
        @(x) isnumeric(x) || islogical(x) || isempty(x));
    parse(p, inputs{:});

    if isempty(p.Results.Upscale)
        upscale = defaultUpscale;
    else
        upscale = p.Results.Upscale;
    end

    if isempty(p.Results.Buffer)
        buffer = defaultBuffer;
    else
        buffer = p.Results.Buffer;
    end

    if isempty(p.Results.Nearby)
        nearby = defaultNearby;
    else
        nearby = logical(p.Results.Nearby);
    end

end
