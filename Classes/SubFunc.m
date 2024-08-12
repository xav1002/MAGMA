classdef SubFunc < handle
    properties
        funcVal = "";
        funcName = "";
        funcSym = "";
        funcLaTeX = "";
        params = {};
        lims = struct('lowerLim',0,'upperLim',0);
    end

    methods
        function subf = SubFunc(funcVal,funcName,funcSym,upperLim,lowerLim)
            subf.funcVal = funcVal;
            subf.funcName = funcName;
            subf.funcSym = funcSym;
            subf.lims.lowerLim = lowerLim;
            subf.lims.upperLim = upperLim;
        end

        % subf: SubFunc class ref, sysVars: string[]
        function initParams(subf,sysVars,defaultParamVals)
            vars = string(findVars(char(subf.funcVal)));
            if strcmp(subf.funcName,'Light Intensity')
                non_params = {'integral','integral2','integral3','l','x','y','z'}';
                if ~isempty(vars)
                    funcParams = setxor(intersect([sysVars(:,2);non_params],vars),vars);
                    funcParams(contains(funcParams,",")) = [];
                else
                    funcParams = [];
                end
            else
                if ~isempty(vars)
                    funcParams = setxor(intersect(sysVars(:,2),vars),vars);
                else
                    funcParams = [];
                end
            end

            prevParams = cell(length(subf.params),2);
            for k=1:1:length(subf.params)
                prevParams{k,1} = convertStringsToChars(subf.params{k}.sym);
                prevParams{k,2} = subf.params{k}.val;
            end
            
            subf.params = {};
            for k=1:1:length(funcParams)
                if ~isempty(prevParams) && any(strcmp(prevParams(:,1),funcParams(k)))
                    paramVal = prevParams{k,2};
                else
                    paramVal = 1;
                end
                subf.updateParams(k,funcParams(k),funcParams(k),paramVal,"",'true',defaultParamVals);
            end
        end

        % subf: SubFunc class ref, name: string, sym: string, val: string, unit: string
        function updateParams(subf,num,name,sym,val,unit,editable,defaultParamVals)
            paramNames = fields(defaultParamVals);
            if strcmp(unit,"mg/L")
                val = val./1000;
            elseif strcmp(unit,"mg/g")
                val = val./1000;
            elseif strcmp(unit,"g/mg")
                val = val.*1000;
            elseif strcmp(unit,"cells/mg")
                val = val.*1000;
            elseif strcmp(unit,"mg/cells")
                val = val./1000;
            end

            if strcmp(sym,"R")
                val = 8.314;
                unit = 'kPa*L/mol*K';
                editable = 'false';
            elseif strcmp(sym,"tau")
                val = 1;
                unit = 's';
                editable = 'false';
            elseif strcmp(sym,"pi")
                val = pi;
                unit = '';
                editable = 'false';
            elseif strcmp(sym,"K_w")
                val = 1E-14;
                unit = 'mol/L';
                editable = 'false';
            elseif any(strcmp(sym,paramNames))
                if contains(sym,'C_p_')
                    unit = 'kJ/kg';
                elseif contains(sym,'rho_')
                    unit = 'g/L';
                elseif contains(sym,'MW_')
                    unit = 'g/mol';
                end
                val = defaultParamVals.(sym);
                editable = 'false';
            end

            newParam = struct( ...
                'num',num, ...
                'func',subf.funcName, ...
                'name',name, ...
                'sym',sym, ...
                'val',val, ...
                'unit',unit, ...
                'editable',editable ...
            );

            addNewParam = true;
            paramIdx = 1;
            for k=1:1:length(subf.params)
                if strcmp(subf.params{k}.sym,sym)
                    addNewParam = false;
                    paramIdx = k;
                    break;
                end
            end

            if addNewParam
                subf.params{end+1} = newParam;
            else
                subf.params{paramIdx} = newParam;
            end
        end

        % subf: SubFunc class ref, varSyms: string[]
        function findParamsInFunc(subf,varSyms,defaultParamVals)
            vars = string(findVars(char(subf.funcVal)));
            paramSyms = setxor(intersect(varSyms(:,2),vars),vars);
            for k=1:1:length(paramSyms)
                % ### FIXME: add functionality to allow user to come back
                % to update envFunc after defining the parameters' values
                funcIdx = strcmp(subf.getSubFuncParamSyms(),paramSyms(k));
                if any(funcIdx)
                    paramVals = subf.getSubFuncParamVals();
                    paramUnits = subf.getSubFuncParamUnits();
                    subf.updateParams(k,paramSyms(k),paramSyms(k),paramVals(funcIdx),paramUnits(funcIdx),'true',defaultParamVals);
                else
                    subf.updateParams(k,paramSyms(k),paramSyms(k),1,"",'true',defaultParamVals);
                end
            end
        end

        % subf: SubFunc class ref, lowerLim: number, upperLim: number
        function updateLims(subf,lowerLim,upperLim)
            subf.lims.lowerLim = lowerLim;
            subf.lims.upperLim = upperLim;
        end

        % subf: SubFunc class ref
        function subfuncName = getSubFuncName(subf)
            subfuncName = subf.funcName;
        end

        % subf: SubFunc class ref
        function subfuncSym = getSubFuncSym(subf)
            subfuncSym = subf.funcSym;
        end

        % subf: SubFunc class ref
        function subfuncVal = getSubFuncVal(subf)
            subfuncVal = subf.funcVal;
        end

        % subf: SubFunc class ref
        function subfuncLaTeX = getSubFuncLaTeX(subf)
            subfuncLaTeX = subf.funcLaTeX;
        end

        % subf: SubFunc class ref
        function paramNames = getSubFuncParamNames(subf)
            paramNames = string.empty(1,0);
            for k=1:1:length(subf.params)
                paramNames(k) = subf.params{k}.name;
            end
        end

        % subf: SubFunc class ref
        function paramNames = getSubFuncParamSyms(subf)
            paramNames = string.empty(1,0);
            for k=1:1:length(subf.params)
                paramNames(k) = subf.params{k}.sym;
            end
        end

        % subf: SubFunc class ref
        function paramVals = getSubFuncParamVals(subf)
            paramVals = zeros(length(subf.params),1);
            for k=1:1:length(paramVals)
                paramVals(k) = subf.params{k}.val;
            end
        end

        % subf: SubFunc class ref
        function paramUnits = getSubFuncParamUnits(subf)
            paramUnits = string.empty(1,0);
            for k=1:1:length(subf.params)
                paramUnits(k) = subf.params{k}.unit;
            end
        end

        % subf: SubFunc class ref, funcName: string
        function setSubFuncName(subf,funcName)
            subf.funcName = funcName;
        end

        % subf: SubFunc class ref, funcVal: string
        function setSubFuncVal(subf,funcVal)
            subf.funcVal = funcVal;
        end

        % subf: SubFunc class ref, funcSym: string
        function setSubFuncSym(subf,funcSym)
            subf.funcSym = funcSym;
        end

        % subf: SubFunc class ref
        function [govFunc,params,reg_param_ct] = compileSubFunc(subf,param_ct,reg_param_ct,regParamList,regression)
            regParamListUpdate = {};
            params = subf.getSubFuncParamVals()';

            govFunc = subf.getSubFuncVal();

            ct = 0;
            reg_param_ct_loc = 0;
            if regression
                for k=1:1:length(subf.params)
                    match = strcmp(regParamList(:,1),char(subf.params{k}.sym));
                    ct = ct + 1;
                    if any(match)
                        if ~strcmp(regParamList{match,2},"")
                            govFunc = replace(govFunc,subf.params{k}.sym,"$("+(regParamList{match,2})+")");
                        else
                            reg_param_ct_loc = reg_param_ct_loc + 1;
                            govFunc = replace(govFunc,subf.params{k}.sym,"$("+(reg_param_ct+reg_param_ct_loc)+")");
                            regParamListUpdate{1,match} = true; %#ok<*AGROW>
                            regParamListUpdate{2,match} = find(match);
                            regParamListUpdate{3,match} = string(reg_param_ct+reg_param_ct_loc);
                        end
                    else
                        govFunc = replace(char(string(govFunc)),subf.params{k}.sym,"#("+(param_ct+ct)+")");
                    end
                end
            else
                syms = subf.getSubFuncParamSyms();
                for l=1:1:length(syms)
                    ct = ct + 1;
                    govFunc = regexprep(char(string(govFunc)),syms(l),"#("+(param_ct+ct)+")");
                end
            end
            reg_param_ct = reg_param_ct + reg_param_ct_loc;
            govFunc = replace(replace(char(string(govFunc)),"#","p"),"$","r");
        end

        % subf: SubFunc class ref, funcLaTeX: string
        function setSubFuncLaTeX(subf,funcLaTeX)
            subf.funcLaTeX = funcLaTeX;
        end

        % subf: SubFunc class ref
        function paramEditable = getSubFuncParamEditable(subf)
            paramEditable = {};
            for k=1:1:length(subf.params)
                paramEditable{k} = subf.params{k}.editable;
            end
        end
    end

    methods (Static)

    end
end