function contrasts = getUniqueContrasts(expInfo)

cc = [];
for ex = 1:length(expInfo)
    cc = [cc expInfo(ex).block.events.contrastValues];
end
contrasts = unique(cc);

end