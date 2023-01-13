 function [reg] = getregion(str)
    s = strfind(str,"CA1");
    if ~(isempty(s)) 
        reg = "CA1";
    end
    s = strfind(str,"CA2");
    if ~(isempty(s)) 
        reg = "CA2";
    end
    s = strfind(str,"CA3");
    if ~(isempty(s)) 
        reg = "CA3";
    end
    s = strfind(str,"DG");
    if ~(isempty(s)) 
        reg = "DG";
    end
 end
    