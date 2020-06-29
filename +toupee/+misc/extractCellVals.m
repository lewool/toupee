function varargout = extractCellVals(c)

for i = 1:numel(c)
    varargout{i} = c{i};
end

end