% func: string
function [vars,idx] = findVars(func)
    func = erase(func,'exp');
    idx = regexp(func,'[^+\-\*/\^()]');

    % ### FIXME: add logic to better detect how to actually deal with
    % underscores
    
    if isempty(idx)
        vars = string(func);
        idx = 1;
        return;
    end

    consecs = {};
    new_arr = [idx(1)];
    for k=1:1:length(idx)
        if k == length(idx)
            consecs{end+1} = new_arr; %#ok<AGROW>
        else
            if idx(k+1) - idx(k) ~= 1
                consecs{end+1} = new_arr; %#ok<AGROW>
                new_arr = [idx(k+1)];
            else
                new_arr(end+1) = idx(k+1); %#ok<AGROW>
                if k == length(idx)-1
                    consecs{end+1} = new_arr; %#ok<AGROW>
                end
            end
        end
    end

    vars = string.empty;
    for k=1:1:length(consecs)
        if any(regexp(char(func(consecs{k})),'[a-zA-Z]')) && ...
            ~strcmp(char(func(consecs{k})),'pi')
            vars(end+1) = string(func(consecs{k})); %#ok<AGROW>
        end
    end
end