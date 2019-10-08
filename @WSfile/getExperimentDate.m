function [experimentDate] = getExperimentDate(obj)

experimentDate = str2num(strcat(obj.file(6:9),obj.file(11:12),obj.file(14:15)));
    
end