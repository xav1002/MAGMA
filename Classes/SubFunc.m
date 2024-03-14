classdef SubFunc < handle
    properties
        funcVal = "";
        funcName = "";
        funcSym = "";
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
        function initParams(subf,sysVars)
            vars = string(symvar(str2sym(subf.funcVal)));
            sysVars(:,2)
            if ~isempty(vars)
                funcParams = setxor(intersect(sysVars(:,2),vars),vars);
            else
                funcParams = [];
            end
            for k=1:1:length(funcParams)
                subf.updateParams(funcParams(k),funcParams(k),1,"");
            end
        end

        % subf: SubFunc class ref, name: string, sym: string, val: string, unit: string
        function updateParams(subf,name,sym,val,unit)
            newParam = struct('name',"",'sym',"",'val',"",'unit',"");
            newParam.name = name;
            newParam.sym = sym;
            newParam.val = val;
            newParam.unit = unit;

            addNewParam = true;
            paramIdx = 1;
            for k=1:1:length(subf.params)
                if subf.params{k}.sym == sym
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
        function findParamsInFunc(subf,varSyms)
            vars = string(symvar(str2sym(subf.funcVal)));
            paramSyms = setxor(intersect(varSyms(:,2),vars),vars);
            for k=1:1:length(paramSyms)
                % ### FIXME: add functionality to allow user to come back
                % to update envFunc after defining the parameters' values
                funcIdx = strcmp(subf.getSubFuncParamSyms(),paramSyms(k));
                if any(funcIdx)
                    paramVals = subf.getSubFuncParamVals();
                    paramUnits = subf.getSubFuncParamUnits();
                    subf.updateParams(paramSyms(k),paramSyms(k),paramVals(funcIdx),paramUnits(funcIdx));
                else
                    subf.updateParams(paramSyms(k),paramSyms(k),1,"");
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
        function subfunc = getSubFuncVal(subf)
            subfunc = subf.funcVal;
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
        function paramEditable = getSubFuncParamEditable(subf)
            paramEditable = {};
            for k=1:1:length(subf.params)
                paramEditable{k} = subf.params{k}.editable; %#ok<AGROW>
            end
        end
    end

    methods (Static)

    end
end