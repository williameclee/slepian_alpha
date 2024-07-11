function XY = joinxy(XYcell)
    % Join the segments (in a cell array) back together, separated by NaNs
    for iCell = 1:length(XYcell) - 1
        XYcell{iCell} = [XYcell{iCell}; [nan, nan]];
    end

    XY = cell2mat(XYcell);
end
