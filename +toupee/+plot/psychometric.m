function f = psychometric(x, y)
end

% Algorithm (pseudocode):
% >> for each sample
%     >> if displacement b/w (sample) : (`tThresh` samples) > `xThresh`
%           >> (sample) : (sample of max displacement) belong to a movement
% >> find all movement starts
% >> for all movements
%     >> merge consecutive, same-direction movements that are separated by 
%     less than `tMinGap` seconds
%     >> find more exact movement start by looking in future (from
%     predefined movement start) for first sample with diff greater than
%     `xOnThresh`
%     >> find more exact end of movement by looking in future (from
%     predefined movement end) for first sample with diff less than 
%     `xOffThresh` 
% >> discard movements that have duration less than `minDur`
% >> get duration of each movement
% >> get total displacement of each movement
% >> get direction of each movement
% >> classify movement ("flinch" if movement duration < `tThresh`)
% >> get movement peak velocity and acceleration