classdef SubFunc < handle
    properties
        funcVal = "";
        funcName = "";
        params = {};
        lims = struct('lowerLim',0,'upperLim',0);
    end

    methods
        function subf = SubFunc(funcVal,funcName,upperLim,lowerLim)
            subf.funcVal = funcVal;
            subf.funcName = funcName;
            subf.lims.lowerLim = lowerLim;
            subf.lims.upperLim = upperLim;
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

        % subf: SubFunc class ref, lowerLim: number, upperLim: number
        function updateLims(subf,lowerLim,upperLim)
            subf.lims.lowerLim = lowerLim;
            subf.lims.upperLim = upperLim;
        end

        % subf: SubFunc class ref
        function subfunc = getSubFuncVal(subf)
            subfunc = subf.funcVal;
        end

        % subf: SubFunc class ref
        function subfuncName = getSubFuncName(subf)
            subfuncName = subf.funcName;
        end

        % subf: SubFunc class ref
        function paramNames = getSubFuncParamSyms(subf)
            paramNames = string.empty;
            for k=1:1:length(paramNames)
                paramNames(k) = subf.params.sym;
            end
        end

        % subf: SubFunc class ref
        function paramVals = getSubFuncParamVals(subf)
            paramVals = zeros(length(subf.params));
            for k=1:1:length(paramVals)
                paramVals(k) = subf.params.val;
            end
        end

        % subf: SubFunc class ref, funcVal: string
        function setSubFuncVal(subf,funcVal)
            subf.funcVal = funcVal;
        end
    end

    methods (Static)

    end
end