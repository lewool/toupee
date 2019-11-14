function [block, Timeline] = loadData(varargin)
% This shit was written by LEW on 30 Jan 2018

%%

if nargin >= 3
	mouseName = varargin{1};
	expDate = varargin{2};
	expNum = varargin{3};
	
	if ischar(expNum)
		expNum = str2double(expNum);
	end 

	expRef = dat.constructExpRef(mouseName,expDate,expNum);
    
elseif nargin == 2
	expRef = varargin{1};
    mouseName = varargin{2};

elseif nargin == 1
	expRef = varargin{1};
end 

currentDataFolder = '\\znas.cortexlab.net\Subjects';
oldDataFolder = '\\zserver.cortexlab.net\Data\Subjects'; 
oldDataFolder2 = '\\zubjects.cortexlab.net\Subjects'; 

try
    blockFilePath = dat.expFilePath(expRef, 'block', 'master');
    load(blockFilePath);
    tlFilePath = dat.expFilePath(expRef, 'timeline', 'master');
    load(tlFilePath);
    
catch
    try
        block = load(strrep(blockFilePath, currentDataFolder, oldDataFolder));
        Timeline = load(strrep(tlFilePath, currentDataFolder, oldDataFolder));
    catch
        try
        block = load(strrep(blockFilePath, currentDataFolder, oldDataFolder2));
        Timeline = load(strrep(tlFilePath, currentDataFolder, oldDataFolder2));
        
        catch %if no TL file
            block = block;
            Timeline = [];
        end
    end

end

end