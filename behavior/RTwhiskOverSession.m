sessionWhisk = mean(eyeData(1).eta.alignedFace{1}(:,91:101,2),2);
sessionRT =  behavioralData(1).wheelMoves.wheelmoves.time' - behavioralData(1).eventTimes(1).daqTime';
subplot(2,1,1); plot(movmean(sessionWhisk,20));
xlabel('Time(s)')
ylabel('Whisking')
subplot(2,1,2); plot(movmean(sessionRT,20));
ylim(0,4)
xlabel('Time(s)')
ylabel('Reaction time') 

%hold on;
%p = polyfit((eyeData(1).eta.eventWindow{1}(:,91:101,2)), sessionWhisk,1);
%yFit = polyval(coefficientsE{irow} , xFit{irow});

