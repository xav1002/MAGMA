classdef Component < handle
    properties
        number = 0; % comps identifying number
        name = ""; % comps name

        transOut = {}; % struct representing mass transport of comps out of system
        initConc = 1; % initial concentration for comps

        funcParams = {}; % cell array of structs that contain the values, names, symbols, and units of parameters,
                        % with the indicies being the parameter numbers.

        sym = "";
        type = "";

        %% funcParams Item Data Structure
        %   funcVal: string
        %   funcName: string
        %   funcCombo: string
        %   funcType: string
        %   reactants?: string (do we need this?)
        %   params: {struct}
        %       locNum: number
        %       gloNum: number
        %       name: string
        %       sym: string
        %       val: number
        %       unit: string
        % lims: {struct - by variable name}
        %       upperVal: number
        %       lowerVal: number

        %% Transport Struct Data Structure
        % trans:
        %   params: number[]
        %   modelType: number that corresponds to the specified transport
        %       model

        rateEq = 1; % overall dX/dt function for comps, all relevant growth and transport terms will be compiled into this function
    end

    methods
        function comp = Component(name, number, initConc, type)
            comp.name = name;
            comp.number = number;
            comp.initConc = initConc;
            if type == "spec"
                comp.sym = "X"+number;
            elseif type == "chem"
                comp.sym = "C"+number;
            end
            comp.type = type;
            
            comp.setDefaultParams(name);
        end

        % comp: comps class ref, name: string
        function comp = setDefaultParams(comp, name)
            
        end

        % comp: comps class ref, prop: string, val: number | string
        function comp = setInitConc(comp, conc)
            comp.initConc = conc;
        end

        % comp: Component class ref
        function initConc = getInitConc(comp)
            initConc = comp.initConc;
        end

        % comp: Component class ref
        function type = getType(comp)
            type = comp.type;
        end

        % comp: comps class ref, funcVal: string, funcName: string,
        % funcCombo: string, funcType: string, paramStr: string{},
        % varNames: {}
        function comp = setModel(comp, funcVal, funcName, funcCombo, funcType, paramStr, varNames, editing, refFuncName)
            if editing
                disp('Func Does Exist')
                idx = comp.getFuncIdx(refFuncName);
                comp.funcParams{idx} = comp.reviseFuncParamObj(funcVal, funcName, funcCombo, funcType, comp.funcParams{idx});
                % creates array of existing param symbols (strings)
                existingParamSyms = cell(1,length(comp.funcParams{idx}.params));
                for k=1:1:length(existingParamSyms)
                    existingParamSyms{k} = comp.funcParams{idx}.params{k}.sym;
                end
                % idxA are the params to be removed, and idxB are the
                % new params to be added
                [~,idxA,idxB] = setxor(existingParamSyms,paramStr);
                for k=1:1:length(idxA)
                    comp.funcParams{idx}.params(idxA(k)) = [];
                end
%                 locNum = comp.funcParams{idx}.params{end}.locNum + 1;
                for k=1:1:length(idxB)
                    comp.funcParams{idx} = comp.createNewParam(1, paramStr{idxB(k)}, paramStr{idxB(k)}, 1, '', comp.funcParams{idx});
%                     locNum = locNum + 1;
                end
%                     paramNums = length(comp.params):1:length(paramNames);
            else
                disp('Func Does not Exist')
                % create new function
                comp.createFuncParamObj(funcVal, funcName, funcCombo, funcType);
                locNum = 1;
                for i=1:1:length(paramStr)
                    comp.funcParams{end} = comp.createNewParam(locNum, paramStr{i}, paramStr{i}, 1, '', comp.funcParams{end});
                    locNum = locNum + 1;
                end
            end
        end

        % comp: comps class ref, func: string, interact: string,
        % chemName: string
        function comp = createFuncParamObj(comp, funcVal, funcName, funcCombo, funcType)
            funcObj = struct( ...
                'funcVal',funcVal, ...
                'funcName',funcName, ...
                'funcCombo',funcCombo, ...
                'funcType', funcType, ...
                'params',struct([]) ...
            );
            comp.funcParams{end+1} = funcObj;
        end
        
        % comp: comps class ref, chemName: string
        function comp = removeFuncParamObj(comp, funcVal, funcName, funcCombo, funcType)
            for k=1:1:length(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcVal,funcVal) && ...
                        strcmp(comp.funcParams{k}.funcName,funcName) && ...
                        strcmp(comp.funcParams{k}.funcCombo,funcCombo) && ...
                        strcmp(comp.funcParams{k}.funcType,funcType)
                    comp.funcParams{k} = [];
                    return;
                end
            end
        end

        % comp: comps class ref, chemName: string, func: string
        function funcExists = checkForDuplicateFunc(comp, funcVal, funcName, funcCombo, funcType)
            funcExists = false;
            for k=1:1:numel(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcVal,funcVal) && ...
                        strcmp(comp.funcParams{k}.funcName,funcName) && ...
                        strcmp(comp.funcParams{k}.funcCombo,funcCombo) && ...
                        strcmp(comp.funcParams{k}.funcType,funcType)
                    funcExists = true;
                    return;
                end
            end
        end

        % comp: comps class ref
        function sym = getSym(comp)
            sym = comp.sym;
        end

        % comp: comps class ref, funcCombo: string, valsOrNames: string,
        % editing: boolean, editedFuncName: string, editedFuncVal: string
        function funcs = getGovFuncs(comp, funcCombo, valsOrNames, editing, refFuncVal, refFuncName)
            % ### FIXME
            funcs = cell(1,2);
            disp(editing+"test")
            segment = 1;
            for k=1:1:length(comp.funcParams)
                if editing
                    if strcmp(comp.funcParams{k}.funcCombo,funcCombo)
                        if valsOrNames == "Names"
                            if strcmp(comp.funcParams{k}.funcName,refFuncName) && strcmp(comp.funcParams{k}.funcVal,refFuncVal)
                                segment = 2;
                            else
                                funcs{segment}{end+1} = comp.funcParams{k}.funcName;
                                disp(k)
                            end
                        elseif valsOrNames == "Values"
                            if strcmp(comp.funcParams{k}.funcName,refFuncName) && strcmp(comp.funcParams{k}.funcVal,refFuncVal)
                                segment = 2;
                            else
                                funcs{segment}{end+1} = comp.funcParams{k}.funcVal;
                            end
                        end
                    end
                else
                    if strcmp(comp.funcParams{k}.funcCombo,funcCombo)
                        if valsOrNames == "Names"
                            funcs{segment}{end+1} = comp.funcParams{k}.funcName;
                        elseif valsOrNames == "Values"
                            funcs{segment}{end+1} = comp.funcParams{k}.funcVal;
                        end
                    end
                end
            end
        end

        % comp: comps class ref, field: string[], val: cell{}
        % returns the indices for comp.funcParams taht correspond to items
        % that match the fields searched
        % make sure that val inputs for text are strings
%         function idx = getFuncByField(comp, fields, val)
%             funcObjs1 = {};
%             for k=1:1:length(comp.funcParams)
%                 for l=1:1:length(fieldnames(comp.funcParams{1}))
%                     funcObjs2 = {};
%                     if comp.funcParams{k}.(field{l}) == val{l}
%                         funcObjs2{l} = comp.funcParams{k}; %#ok<AGROW> 
%                     end
%                     funcObjs1{k} = funcObjs2; %#ok<AGROW> 
%                 end
%             end
% 
%             funcObjs = intersect(funcObjs1{1},funcObjs1{2});
%             for k=3:1:length(funcObjs1)
%                 funcObjs = intersect(funcObjs,funcObjs1{k});
%             end

%             comboFuncParam = comp.funcParams{:};
%             idxArr = zeros(length(fields),length(comp.funcParams))
%             disp(comboFuncParam)
%             disp([comboFuncParam.(fields(1))])
%             disp(val{1})
%             for k=1:1:length(fields)
%                 if isa(val{k},'char')
%                     val{k} = convertCharsToStrings(val{k});
%                 end
%                 disp([comboFuncParam.(fields(k))])
%                 idxArr(k,:) = ([comboFuncParam.(fields(k))] == val{k});
%             end
%             if length(fields) < 2
%                 idx = idxArr;
%             end
%             for k=1:1:(length(fields)-1)
%                 idx = intersect(idxArr(k,:),idxArr(k+1,:));
%             end
%         end

        % comp: comps class ref, paramSym: string, newVal: number,
        % newUnit: string, newParamName: string, funcName: string
        function comp = updateParam(comp, paramSym, newVal, newUnit, newParamName, funcName)
            funcObjIdx = comp.getFuncIdx(funcName);
            comp.funcParams{funcObjIdx} = comp.reviseParam(paramSym, newVal, newUnit, newParamName, comp.funcParams{funcObjIdx});
        end

        % comp: comps class ref, funcName: string,
        function idx = getFuncIdx(comp, funcName)
            % ### FIXME: if functions have the same value, then there is an
            % error because the first matching function is pulled
            for k=1:1:length(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcName,funcName)
                    idx = k;
                    return;
                end
            end
            disp("No Function Found, Component.m, line 223");
        end

        % comp: comps class ref
        % returns parameter info ready to display to ODESys (params struct)
        function params = getGrthParams(comp)
            params = cell([comp.getParamCount(),6]);
            ct = 1;
            for k=1:1:length(comp.funcParams)
                for j=1:1:length(comp.funcParams{k}.params)
                    params{ct,1} = num2str(ct);
                    params{ct,2} = comp.funcParams{k}.funcName;
                    params{ct,3} = comp.funcParams{k}.params{j}.sym;
                    params{ct,4} = comp.funcParams{k}.params{j}.name;
                    params{ct,5} = comp.funcParams{k}.params{j}.val;
                    params{ct,6} = comp.funcParams{k}.params{j}.unit;
                    ct = ct + 1;
                end
            end
        end

        % comp: comps class ref
        function paramCt = getParamCount(comp)
            paramCt = 0;
            for k=1:1:length(comp.funcParams)
                paramCt = paramCt + length(comp.funcParams{k}.params);
            end
        end

        % comp: comps class ref
        function removeModel(comp,funcVal,funcName,funcCombo,funcType)
            for k=1:1:length(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcVal,funcVal) && ...
                        strcmp(comp.funcParams{k}.funcName,funcName) && ...
                        strcmp(comp.funcParams{k}.funcCombo,funcCombo) && ...
                        strcmp(comp.funcParams{k}.funcType,funcType)
                    comp.funcParams(k) = [];
                    return;
                end
            end
        end

        % comp: comps class ref, funcCombo: string
        function funcNames = getFuncNames(comp,funcCombo)
            funcNames = {};
            ct = 1;
            for k=1:1:length(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcCombo, funcCombo)
                    funcNames{ct} = comp.funcParams{k}.funcName; %#ok<AGROW> 
                    ct = ct + 1;
                end
            end
        end

        % comp: Component class ref, funcName: string
        function funcVal = getFuncValByName(comp,funcName)
            % ### FIXME
            for k=1:1:length(comp.funcParams)
%                 comp.funcParams{k}
                if strcmp(comp.funcParams{k}.funcName, funcName)
                    funcVal = comp.funcParams{k}.funcVal;
                    return;
                end
            end
        end

        % comp: Component class ref
        function funcNum = getFuncCount(comp)
            funcNum = length(comp.funcParams);
        end

        % comp: Component class ref, paramCt: numbers, regParamList: {}, regression: boolean
        function [govFunc,params,reg_param_ct] = compileGovFunc(comp,param_ct,reg_param_ct,regParamList,regression)
            % Need to make sure that each individual govFunc component is
            % isolated (so (C1+C2)*(C3+C4) ~= C1+C2*C3_C4)
            govFunc = "";
            params = [];
            for k=1:1:length(comp.funcParams)
                if comp.funcParams{k}.funcCombo == "Multiplicative"
                    if k == 1
                        govFunc = govFunc + comp.funcParams{k}.funcVal;
                    else
                        govFunc = govFunc + "*" + comp.funcParams{k}.funcVal;
                    end
                end
                % compiling parameters array
                for l=1:1:length(comp.funcParams{k}.params)
                    params = [params,comp.funcParams{k}.params{l}.val]; %#ok<AGROW>
                end
            end

            if govFunc == ""
                operator = "";
            elseif govFunc ~= ""
                operator = "+";
            end
            for k=1:1:length(comp.funcParams)
                if comp.funcParams{k}.funcCombo == "Additive"
                    govFunc = govFunc + operator + comp.funcParams{k}.funcVal;
                end
            end

            % replacing param syms in govFunc with p({gloNum})
            ct = 1;
            if regression
                for k=1:1:length(comp.funcParams{k})
                    for l=1:1:length(comp.funcParams{k}.params)
                        match = strcmp(cellstr(regParamList),char(comp.funcParams{k}.params{l}.sym));
                        if any(match)
                            govFunc = split(govFunc,"#");
                            govFunc = govFunc(strlength(govFunc) > 1);
                            govFunc = regexprep(govFunc,comp.funcParams{k}.params{l}.sym,"#r#("+(reg_param_ct)+")");
                            reg_param_ct = reg_param_ct + 1;
                        else
                            govFunc = split(govFunc,"#");
                            govFunc = govFunc(strlength(govFunc) > 1);
                            govFunc = regexprep(govFunc,comp.funcParams{k}.params{l}.sym,"#p#("+(param_ct+ct)+")");
                            ct = ct + 1;
                        end
                    end
                end
            else
                for k=1:1:length(comp.funcParams{k})
                    for l=1:1:length(comp.funcParams{k}.params)
                        govFunc = split(govFunc,"#");
                        govFunc = govFunc(strlength(govFunc) > 1);
                        govFunc = regexprep(govFunc,comp.funcParams{k}.params{l}.sym,"p("+(param_ct+ct)+")");
                        ct = ct + 1;
                    end
                end
            end
            govFunc = strrep(govFunc,"#","");
        end
    end

    methods (Static)
        % num: number, name: string, sym: string,
        % val: number | string, unit: string, funcObj: struct
        function funcObj = createNewParam(locNum, name, sym, val, unit, funcObj)
            % won't maintain global number until very end, when all
            % comps and chemicals are set
            param = struct( ...
                'locNum',locNum, ...
                'gloNum',1, ...
                'name',name, ...
                'sym',sym, ...
                'val',val, ...
                'unit',unit ...
            );
            funcObj.params{end+1} = param;
        end

        % paramSym: string, newVal: number, funcObj: struct
        function funcObj = reviseParam(paramSym, newVal, newUnit, newParamName, funcObj)
            paramSyms = cell(1,numel(funcObj.params));
            for k=1:1:numel(funcObj.params)
                paramSyms{k} = funcObj.params{k}.sym;
            end
            idx = strcmp(paramSyms,paramSym);
            funcObj.params{idx}.name = newParamName;
            funcObj.params{idx}.val = newVal;
            funcObj.params{idx}.unit = newUnit;
        end

        % paramName: string, funcObj: struct
        function funcObj = removeParam(paramName, funcObj)
            for k=1:1:length(funcObj.params)
                if strcmp(funcObj.params{k}.name, paramName)
                    % removes parameter from funcObj
                    funcObj.params{k} = [];
                    return;
                end
            end
        end

        % funcVal: string, funcName: string, funcCombo: string, funcType: string, funcObj: struct
        function newFuncObj = reviseFuncParamObj(funcVal, funcName, funcCombo, funcType, funcObj)
            funcObj.funcVal = funcVal;
            funcObj.funcName = funcName;
            funcObj.funcCombo = funcCombo;
            funcObj.funcType = funcType;
            newFuncObj = funcObj;
        end
    end
end