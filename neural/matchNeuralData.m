function matched = matchNeuralData(neuralData, expInfo, behavioralData)

for a = 1:length(neuralData(1).eta.alignedResps)
    for s = 1:length(neuralData)
        hc{s,a} = neuralData(s).eta.alignedResps{a};
    end
end

for a = 1:length(neuralData(1).eta.alignedResps)
    matched.neuralData.eta.alignedResps{a} = cat(1,hc{:,a});
end

matched.neuralData.eta.events = neuralData(1).eta.events;
matched.neuralData.eta.eventWindow = neuralData(1).eta.eventWindow;

matched.neuralData = getSignificantActivity(expInfo, behavioralData, matched.neuralData,1);
