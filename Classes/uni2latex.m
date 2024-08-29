function latex_out = uni2latex(eqn)
    % capturing first operator (usually negative)
    if length(eqn) == 1
        latex_out = eqn;
        return;
    end

    num_left_paren = regexp(eqn,'\(');
    num_right_paren = regexp(eqn,'\)');
    if length(num_left_paren) ~= length(num_right_paren), return; end
    if regexp(eqn(end),'[+\-\*/^=]','once'), return; end
    % if strcmp(eqn(1),'-') && ~strcmp(eqn(2),'(')
    %     leadOperator = '-';
    % else
    %     leadOperator = '';
    % end

    % correcting exp()
    while ~isempty(regexp(eqn,'exp(','once'))
        % use the parenthesis logic here
        exp_idx = regexp(eqn,'exp(','once');
        operator_idx = regexp(eqn,'[+\-\*/^=\(\)]');
        next_operator_idx = operator_idx(find((exp_idx+4 < operator_idx),1));
        if isempty(next_operator_idx)
            eqn = [eqn(1:exp_idx-1),'e^(',eqn(exp_idx+4:end),')'];
        else
            num_left = 1;
            num_right = 0;
            right_eqn = eqn(exp_idx+4:end);
            end_paren_idx = exp_idx+4;
            while num_left ~= num_right
                next_paren = regexp(right_eqn,'[()]','once');
                if strcmp(right_eqn(next_paren),')')
                    num_right = num_right + 1;
                elseif strcmp(right_eqn(next_paren),'(')
                    num_left = num_left + 1;
                end
                right_eqn = right_eqn(next_paren+1:end);
                end_paren_idx = end_paren_idx+next_paren;
            end

            next_operator_idx = operator_idx(find((end_paren_idx-1 <= operator_idx),1));
            if next_operator_idx == operator_idx(end)
                eqn = [eqn(1:exp_idx-1),'e^{(',eqn(exp_idx+4:end),'}'];
            else
                eqn = [eqn(1:exp_idx-1),'e^({',eqn(exp_idx+4:next_operator_idx-1),'}',eqn(next_operator_idx:end)];
            end
        end
    end
    
    % correcting subscripts
    while ~isempty(regexp(eqn,'_[^{]','once'))
        underscore_idx = regexp(eqn,'_[^{]','once');
        operator_idx = regexp(eqn,'[+\-\*/^=()]');
        next_operator_idx = operator_idx(find((underscore_idx < operator_idx),1));

        if isempty(next_operator_idx)
            eqn = [eqn(1:underscore_idx),'{',eqn(underscore_idx+1:end),'}'];
        else
            if strcmp(eqn(underscore_idx+1),'(')
                num_left = 0;
                num_right = 1;
                right_eqn = eqn(underscore_idx+2:end);
                while num_left ~= num_right
                    next_paren = regexp(right_eqn,'[()]','once');
                    if strcmp(right_eqn(next_paren),')')
                        num_right = num_right + 1;
                    elseif strcmp(right_eqn(next_paren),'(')
                        num_left = num_left + 1;
                    end
                    disp([num_left,num_right])
                end

                end_paren_idx = next_paren;
                next_operator_idx = operator_idx(find((end_paren_idx < operator_idx),1))+1;
            end
            eqn = [eqn(1:underscore_idx),'{',eqn(underscore_idx+1:next_operator_idx-1),'}',eqn(next_operator_idx:end)];
        end
    end

    % correcting fractions
    while ~isempty(regexp(eqn,'/','once'))
        div_idx = regexp(eqn,'/','once');
        operator_idx = regexp(eqn,'[+\-\*/^=\(\)]');
        % ### EDGE: causes inaccuracies with parenthesis if divide is
        % wrapped in parenthesis

        next_operator_idx = operator_idx(find((div_idx < operator_idx),1));

        if isempty(next_operator_idx)
            if strcmp(eqn(div_idx-1),')')
                num_left = 0;
                num_right = 1;
                left_eqn = reverse(eqn(1:div_idx-2));
                left_eqn_new = left_eqn;
                while num_left ~= num_right
                    next_paren = regexp(left_eqn_new,'[\(\)]','once');
                    if strcmp(left_eqn_new(next_paren),')')
                        num_right = num_right + 1;
                    elseif strcmp(left_eqn_new(next_paren),'(')
                        num_left = num_left + 1;
                    end
                    left_eqn_new = left_eqn_new(next_paren+1:end);
                end

                end_paren_idx = next_paren;
                if end_paren_idx == 0
                    eqn = ['\frac{',eqn(1:operator_idx(end)-1),'}{',eqn(operator_idx(end)+1:end),'}'];
                else
                    prev_operator_idx = operator_idx(find((end_paren_idx < operator_idx),1))-1;
                    eqn = [eqn(1:prev_operator_idx+1),'\frac{',eqn(prev_operator_idx+2:operator_idx(end)-1),'}{',eqn(operator_idx(end)+1:end),'}'];
                end
            elseif length(operator_idx) == 1
                eqn = ['\frac{',eqn(1:operator_idx(end)-1),'}{',eqn(operator_idx(end)+1:end),'}'];
            else
                eqn = [eqn(1:operator_idx(end-1)),'\frac{',eqn(operator_idx(end-1)+1:operator_idx(end)-1),'}{',eqn(operator_idx(end)+1:end),'}'];
            end
        else
            if div_idx == operator_idx(1)
                prev_operator_idx = 1;
            else
                prev_operator_idx = operator_idx(find((div_idx < operator_idx),1)-2);
            end

            if strcmp(eqn(div_idx-1),')')
                num_left = 0;
                num_right = 1;
                left_eqn = reverse(eqn(1:div_idx-2));
                left_eqn_2 = reverse(eqn(1:div_idx-2));
                while num_left ~= num_right
                    next_paren = regexp(left_eqn_2,'[()]','once');
                    if strcmp(left_eqn(next_paren),')')
                        num_right = num_right + 1;
                        left_eqn_2 = regexprep(left_eqn_2,')',' ','once');
                    elseif strcmp(left_eqn(next_paren),'(')
                        num_left = num_left + 1;
                        left_eqn_2 = regexprep(left_eqn_2,'(',' ','once');
                    end
                end

                end_paren_idx = length(left_eqn) + 1 - next_paren;
                prev_operator_idx = end_paren_idx - 1;
            end

            end_paren_idx = 1;
            is_right_paren = 1;
            if strcmp(eqn(div_idx+1),'(')
                num_left = 1;
                num_right = 0;
                right_eqn = eqn(div_idx+2:end);
                right_eqn_2 = eqn(div_idx+2:end);
                while num_left ~= num_right
                    next_paren = regexp(right_eqn_2,'[()]','once');
                    if strcmp(right_eqn(next_paren),')')
                        num_right = num_right + 1;
                        right_eqn_2 = regexprep(right_eqn_2,')',' ','once');
                    elseif strcmp(right_eqn(next_paren),'(')
                        num_left = num_left + 1;
                        right_eqn_2 = regexprep(right_eqn_2,'(',' ','once');
                    end
                end

                end_paren_idx = next_paren + div_idx;
                next_operator_idx = operator_idx(find((end_paren_idx < operator_idx),1));
                is_right_paren = 0;
            end

            if div_idx == operator_idx(1)
                eqn = ['\frac{',eqn(1:div_idx-1),'}{', ...
                    eqn(div_idx+1:next_operator_idx-is_right_paren),'}',eqn(next_operator_idx+(1-is_right_paren):end)];
            elseif end_paren_idx == operator_idx(end)
                eqn = [eqn(1:prev_operator_idx),'\frac{',eqn(prev_operator_idx+1:div_idx-1),'}{', ...
                eqn(div_idx+1:end),'}'];
            else
                eqn = [eqn(1:prev_operator_idx),'\frac{',eqn(prev_operator_idx+1:div_idx-1),'}{', ...
                eqn(div_idx+1:next_operator_idx-is_right_paren),'}',eqn(next_operator_idx+(1-is_right_paren):end)];
            end
        end
    end

    % if the operator directly next to the underscore, divide, or
    % exponent is ( or ), then need extra logic to recursively find the
    % closing parenthesis

    eqn = replace(eqn,'(','\left(');
    eqn = replace(eqn,')','\right)');
    % ### FIXME: how to not take greek characters that are part of words?
    eqn = regexprep(eqn,'tau','\\tau');
    eqn = regexprep(eqn,'mu','\\mu');
    eqn = regexprep(eqn,'chi','\\chi');
    eqn = regexprep(eqn,'beta','\\beta');
    eqn = regexprep(eqn,'alpha','\\alpha');
    % eqn = [leadOperator,eqn];
    latex_out = replace(eqn,'*','\,');
end