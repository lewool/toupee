function contrasts = getAllContrasts(expInfo)

cc = [];
for ex = 1:length(expInfo)
    cc = [cc expInfo(ex).block.events.contrastValues];
end
contrasts = cc;

end