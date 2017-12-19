function p = read_cfg(cfgFn)
    fid = fopen(cfgFn);
    if fid == -1
        disp('invalid file');
        p = 0;
        return
    end
    tline = fgetl(fid);
    pair = cell(50,1);
    nl = 1;
    while tline(1) ~= -1
        if tline(1) ~= '#'
            element = strtrim(strsplit(tline,'='));
            pair{nl} = element{1};
            noComments = strtrim(strsplit(element{2},' '));
            element(2) = noComments(1);
            match = regexp(element{2},'(^\-?\.?\d*$|^\-?\d*\.?\d*$)','once');
            if ~isempty(match)
                pair{nl+1} = str2double(element{2});
            else
                logicals = regexp(element{2},'(^(T|t)rue$|^(F|f)alse$)','once');
                if logicals
                    if isequal(element{2},'True') || isequal(element{2},'true') || isequal(element{2},'TRUE')
                        pair{nl+1} = true;
                    else
                        pair{nl+1} = false;
                    end
                else
                    pair{nl+1} = element{2};
                end
            end
            multiple = false;
            for i=nl-2:-2:1
                if isequal(pair{nl}, pair{i})
                    pair{i+1} = [pair{i+1}, pair{nl+1}];
                    multiple = true;
                    break
                end
            end
            if ~multiple
                nl = nl +2;
            end
        end
        tline = fgetl(fid);
    end
    pair = pair(1:nl-1);
    p = struct(pair{:});
end
