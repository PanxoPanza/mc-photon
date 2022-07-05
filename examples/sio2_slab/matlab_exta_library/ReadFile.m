function Data  = ReadFile(FileName,varargin)
warning off backtrace

% process instructions
paramMatch = {'Verbose','on'};
rule = parseInput(paramMatch, varargin{:});


if strcmp(rule.Verbose,'on')
    disp(['Opening: ' FileName '...'])
end
SETfid = fopen(FileName);
if SETfid ~=-1
    Ncol_max = 0;
    % get the number of lines in the file
    nLines = 0;
    odata = ftell(SETfid); % get the first line of data
    while (~feof(SETfid))
        tmp = fgetl(SETfid);
        tmp = strtrim(tmp);
        if ~isempty(tmp) ...            % exclude lines with comments
                && ~strcmp(tmp(1),'#') ...
                && ~strcmp(tmp(1),'%') ...
                && ~strcmp(tmp(1),'//')
            Ncol = length(strsplit(tmp,{'\t', ' ', ',','"'})); % number of columns
            if Ncol > Ncol_max
                Ncol_max = Ncol;
            end
            nLines = nLines + 1;
        else
            odata = ftell(SETfid);
        end
    end
    fseek(SETfid,odata,'bof'); % return pointer to the first line of data
    
    % extract data
    Data = zeros(nLines,Ncol_max);
    for idata = 1:nLines
        tmp = fgetl(SETfid);
        tmp = strtrim(tmp);
        output  = strsplit(tmp,{'\t', ' ', ',','"'});
        icol = 0;
        for icol_data = 1:length(output)
            if ~isempty(output{icol_data})
                icol = icol + 1;
                Data(idata,icol) = str2double(output{icol_data});
            end
        end
    end
    if strcmp(rule.Verbose,'on')
        disp('... File successfully open')
    end
    fclose(SETfid);
    return
else
    warning(['File: ' FileName ' not found!'])
    Data = 0;
    return
end
end