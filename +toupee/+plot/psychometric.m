function f = psychometric(x, y)
end

% Algorithm (pseudocode):
% >> for each sample
%     >> if displacement b/w (sample) : (`tThresh` samples) > `xThresh`
%           >> (sample) : (sample of max displacement) belong to a movement
% >> find all movement starts
% >> merge consecutive, same-direction movements that are separated by 
%     less than `tMinGap` seconds
% >> for all movements
%     >> get direction of each movement
%     >> find more exact movement start by looking in future (from
%     predefined movement start) for first sample with diff greater than
%     `xOnThresh`
%     >> find more exact end of movement by looking in future (from
%     predefined movement end) for first sample with diff less than 
%     `xOffThresh` 
% >> get duration of each movement
% >> discard movements that have duration less than `minDur`
% >> get total displacement of each movement
% >> classify movement ("flinch" if movement duration < `tThresh`)
% >> get movement peak velocity and acceleration

for iM = 1:15
    figure
    plot(startS(mergeeMove(iM):(mergeeMove(iM) + 20)), x(startS(mergeeMove(iM)):(startS(mergeeMove(iM)) + 20)))
    hold on
    plot(startS(mergerMove(iM):(mergerMove(iM) + 20)), x(startS(mergerMove(iM)):(startS(mergerMove(iM)) + 20)))
    title(sprintf('%i: Concord %i, Diff %i', iM, concordMovesMask(iM), startS(mergerMove(iM)) - startS(mergeeMove(iM))))
    legend('mergee', 'merger')
end

%     % If there is a change in direction, just keep the movement in the
%     % first direction
%     % If the samples pass thresh in both directions, just keep the move in
%     % the first direction, the second will be caught later.
%     % direction of displacement for current samples
%     curDisDir = sign(diff(curDis(curMoveMask)));
%     if numel(unique(curDisDir)) == 3
%         iFirstMoveStop = find(curDisDir == -curDisDir(1), 1, 'first');
%         curMoveIdxs = find(curMoveMask);
%         curMoveMask(curMoveIdxs(iFirstMoveStop:end)) = 0;
%     end
%     % If a movement that passes thresh has position change in the opposite
%     % direction beforehand, don't count it; it will be caught later.
%     iEndMove = find(curMoveMask, 1, 'last') - 1;
%     % direction of position change for current move
%     curMoveDir = sign(diff(x(iS:(iS + iEndMove))));
%     if numel(find(curMoveDir == -1)) > 2 ...
%        && numel(find(curMoveDir == 1)) > 2 
%         continue
%     % Else, keep the movement and mark these samples.
%     else
%         moveDir(iS:(iS + iEndMove)) = true;
%     end