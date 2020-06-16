function [expRef, expLog] = constructExpRef(mouseName,expDate,expNum)
    expRef = strcat(expDate,'_',num2str(expNum),'_',mouseName);
    expLog = {{mouseName};{expDate};{num2str(expNum)}};
end