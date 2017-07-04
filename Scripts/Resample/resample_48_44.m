function [downsampleFile] = resample_48_44(file)
        upsampleFile = interp(file,147);
        downsampleFile = downsample(upsampleFile,160);
    end