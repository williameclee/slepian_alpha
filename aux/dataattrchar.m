function dataFileAttr = dataattrchar(varargin)
    p = inputParser;
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 0);
    addOptional(p, 'Inclang', 90, ...
        @(x) isnumeric(x) && length(x) <= 2);
    addOptional(p, 'EQMask', 0, ...
        @(x) isnumeric(x));
    addOptional(p, 'MoreBuffer', []);
    parse(p, varargin{:});
    upscale = p.Results.Upscale;
    inclang = p.Results.Inclang;
    buf = p.Results.Buffer;
    eqMask = p.Results.EQMask;
    moreBuf = p.Results.MoreBuffer;

    %% Find the data file
    dataFileAttr = cell(1, 3);

    if ~(upscale == 0 || upscale == 1 || isempty(upscale))
        dataFileAttr{1} = num2str(upscale);
    end

    if ~(buf == 0) || ~isempty(moreBuf)
        dataFileAttr{2} = num2str(buf);

        if ~isempty(moreBuf)

            for i = 1:length(moreBuf)
                moreBuf{i} = char(join(string([moreBuf{i}]), ''));
            end

            moreBuf = char(join(moreBuf, '_'));

            dataFileAttr{2} = [dataFileAttr{2}, '_', moreBuf];
        end

    end

    if isnan(inclang) || isempty(inclang)
    elseif isscalar(inclang)

        if ~(inclang == 90 | isempty(inclang))
            dataFileAttr{3} = num2str(inclang);
        end

    else
        inclangSign = [];

        if inclang(1) < 0
            inclangSign(1) = 's';
        elseif inclang(1) > 0
            inclangSign(1) = 'n';
        end

        if inclang(2) < 0
            inclangSign(2) = 's';
        elseif inclang(2) > 0
            inclangSign(2) = 'n';
        end

        dataFileAttr{3} = ...
            [num2str(abs(inclang(1))), inclangSign(1), '_', ...
             num2str(abs(inclang(2))), inclangSign(2)];

    end

    % Find the last filled attribute
    iLastUsedAttr = find(~cellfun(@isempty, dataFileAttr), 1, 'last');

    if ~isempty(iLastUsedAttr)
        dataFileAttr = dataFileAttr(1:iLastUsedAttr);
        % Fill in 0
        emptyBeforeLast = 1:(iLastUsedAttr - 1);
        emptyCellsIndex = cellfun(@isempty, dataFileAttr(emptyBeforeLast));
        dataFileAttr(emptyBeforeLast(emptyCellsIndex)) = {'0'};
        % Add a hyphen before each attribute
        if length(dataFileAttr) > 1
            dataFileAttr(2:end) = cellfun(@(x) ['-', x], dataFileAttr(2:end), ...
                'UniformOutput', false);
        end

    end

    if all(cellfun(@isempty, dataFileAttr))
        dataFileAttr = '';
    else
        dataFileAttr = char(join(dataFileAttr, ''));
    end

    if eqMask ~= 0

        if isempty(dataFileAttr)
            dataFileAttr = ['EQ', num2str(eqMask)];
        else
            dataFileAttr = [dataFileAttr, '-EQ', num2str(eqMask)];
        end

    end

end
