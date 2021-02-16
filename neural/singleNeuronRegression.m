function singleNeuronRegression(X,Y)
% X is an array of neurons x trials, like one of the output arrays from
% getEpochResps.m where there is one value per neuron per trial.
% Y is a vector of outcomes/IDs, such as contrasts or choices.

for iX = 1:length(X)
    