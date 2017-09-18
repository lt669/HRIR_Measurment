function [] = ITD_Lookup_Table_Generation(projectName,subjectName,fs)
    disp('');    
    disp('---Running IDT_Lookup_Table_Generation---');
    
    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    
    % ---LOAD FILE PATHS--- %

    hrirPath = strcat('Audio/Projects/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder);
    hrirDirectory = dir(sprintf('%s/*.wav',hrirPath));

    % --------------------- %

    for k =1:length(hrirDirectory)

        fileName = hrirDirectory(k).name;

        % Scan and store angles
        [azi_angle,ele_angle] = find_azi_ele(fileName);
        
        itd_out(k,1) = azi_angle;
        itd_out(k,2) = ele_angle;
        
        % Load audio file 
        hrir = audioread(sprintf('%s/%s',hrirPath,fileName));

        % Find ITD
        maxLag = round(0.0011*fs);

        left = hrir(:,1);
        right = hrir(:,2);

        itd = finddelay(left,right,maxLag);
        itd_out(k,3) = (itd/fs)*10^3;

        itd_out_filepath = strcat('Audio/Projects/',projectName,'/SOFAFiles/lookUpTables');
        disp(strcat('ITD: ',itd_out_filepath));
        if 0==exist(itd_out_filepath,'file')
            disp('Creating File...');
            disp(itd_out_filepath);
            mkdir(itd_out_filepath);
        end
        % Write to text file
        idt_out_name = strcat(itd_out_filepath,'/',subjectName,'_ITD_Lookup_Table.txt');

        fid = fopen(idt_out_name,'w');
        fprintf(fid,'%0i \t %0i \t %0f\r\n',itd_out');
        fclose(fid);
        
    end
    
%     plotOut(1) = itd_out(3,3);
%     plotOut(2) = itd_out(43,3);
%     plotOut(3) = itd_out(49,3);
%     plotOut(4) = itd_out(10,3);
%     plotOut(5) = itd_out(16,3);
%     plotOut(6) = itd_out(24,3);
%     plotOut(7) = itd_out(30,3);
%     plotOut(8) = itd_out(36,3);
%     plotOut(8) = itd_out(3,3);
%         disp('Plotting ITD...');
        % Plot ITD
       % plot(plotOut);
    
    function [azi,ele] = find_azi_ele(name)
        [token,remain] = strtok(name,'_');
        [azi,remain] = strtok(remain,'_');
        [ele,remain] = strtok(remain,'_ele');
        
        azi = str2num(azi);
        ele = str2num(ele);
        
        
        
    end
disp('');    
disp('---Exiting IDT_Lookup_Table_Generation---');
end