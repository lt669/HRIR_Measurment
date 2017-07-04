%{
    Function to find and output the azimuth and elevation file name without
    any other extensions.
%}

function [outputName] = findName(fullPath)
    
[pathstr,inputName,ext] = fileparts(fullPath);

filenamestr = char(inputName); % Get current filename

IndexAzi = strfind(inputName, 'azi_'); % Find the text 'azi_'
Azout = sscanf(filenamestr(1,IndexAzi + ...
    length('azi_'):end), '%g', 1); % Get azimuth value

IndexEle = strfind(inputName, 'ele_'); % Find the text 'ele_'
Elout = sscanf(filenamestr(1,IndexEle + ...
    length('ele_'):end), '%g', 1); % Get elevation value

Azout = int2str(Azout);
Elout = int2str(Elout);

outputName = cellstr(strcat('azi_',Azout,'_ele_',Elout));
disp(sprintf('findName: outputName = %s',char(outputName)));

end