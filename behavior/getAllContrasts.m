function contrasts = getAllContrasts(expInfo)

cc = [];
for ex = 1:length(expInfo)
    nt = length(expInfo(ex).block.events.endTrialValues);
    cc = [cc expInfo(ex).block.events.contrastValues(1:nt)];
end
contrasts = cc;

end