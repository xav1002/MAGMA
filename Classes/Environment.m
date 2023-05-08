classdef Environment < handle
    properties
        name = '';

        subfuncs = {};
    end

    methods
        function env = Environment(name,modelVarSyms)
            env.name = name;

            env.setDefault(name,modelVarSyms);
        end

        % env: Environment class ref, val: number
        function setModelTime(env, val)
            env.subfuncs{1}.setSubFuncVal(val);
        end

        % env: Environment class ref, funcStr: string
        function setLightFunc(env, funcStr)
            env.subfuncs{2}.setSubFuncVal(funcStr);
        end

        % env: Environment class ref, funcStr: string
        function setTempFunc(env, funcStr)
            env.subfuncs{3}.setSubFuncVal(funcStr);
        end

        % env: Environment class ref, val: number
        function setCulVol(env, val)
            env.subfuncs{4}.setSubFuncVal(val);
        end

        % env: Environment class ref, val: number
        function setCulSA(env, val)
            env.subfuncs{5}.setSubFuncVal(val);
        end

        % env: Environment class ref, envFuncName: string, paramSym: string,
        % newVal: number, newUnit: string, newParamName: string
        function updateEnvFuncParams(env,envFuncName,paramSym,newVal,newUnit,newParamName)
            for k=1:1:length(env.subfuncs)
                if strcmp(env.subfuncs{k}.getSubFuncName(),envFuncName)
                    env.subfuncs{k}.updateParams(newParamName,paramSym,newVal,newUnit);
                    break;
                end
            end
        end

        % env: Environment class ref, envFuncName: string
        function params = getGrthParamsByEnvFuncName(env,envFuncName)
            for k=1:1:length(env.subfuncs)
                if strcmp(env.subfuncs{k}.getSubFuncName(),envFuncName)
                    paramNames = env.subfuncs{k}.getSubFuncParamNames();
                    paramSyms = env.subfuncs{k}.getSubFuncParamSyms();
                    paramVals = env.subfuncs{k}.getSubFuncParamVals();
                    paramUnits = env.subfuncs{k}.getSubFuncParamUnits();
                    paramNums = 1:1:length(paramNames);
                    params = cell(length(paramNums),6);
                    for l=paramNums
                        params{l,1} = char(string(paramNums(l)));
                        params{l,2} = char(envFuncName);
                        params{l,3} = char(paramSyms(l));
                        params{l,4} = char(paramNames(l));
                        params{l,5} = paramVals(l);
                        params{l,6} = char(paramUnits(l));
                    end
                end
            end
        end

        % env: Environment class ref, name: string, modelVarSyms: string[]
        function setDefault(env, name, modelVarSyms)
            envDef = EnvDefaults();
            if ~any(envDef.Names == name)
                name = "Custom_Env";
            end
            envParamVals = {};
            envParamVals{1} = envDef.(name).modelTime;
            envParamVals{2} = envDef.(name).lightFunc;
            envParamVals{3} = envDef.(name).tempFunc;
            envParamVals{4} = envDef.(name).culVol;
            envParamVals{5} = envDef.(name).culSA;
            envParamNames = env.getParamNames();
            for k=1:1:size(envParamNames,1)
                env.subfuncs{k} = SubFunc(envParamVals{k},envParamNames{k,1},envParamNames{k,2},0,0);
                env.subfuncs{k}.initParams(modelVarSyms);
%                 env.subfuncs{k}.updateParams(envParamNames{k,2},envParamNames{k,2},1,"");
            end
        end

        % env: Environment class ref
        function paramVals = getParamVals(env)
            paramVals = cell(size(env.subfuncs'));
            for k=1:1:length(paramVals)
                func = env.subfuncs{k}.getSubFuncVal();
                paramVals{k} = char(func);
            end
        end

        % env: Environment Class ref, varSyms: {}
        function funcHandles = getParamFuncHandles(env,varSyms)
            funcHandles = cell(size(env.subfuncs'));
            for k=1:1:length(funcHandles)
                func = env.subfuncs{k}.getSubFuncVal();
                env.subfuncs{k}.findParamsInFunc(varSyms);
                func = env.evalParams(func);
                funcHandles{k} = char(func);
            end
        end

        % env: Environment class ref, func: string
        function func = evalParams(env,func)
            for k=1:1:length(env.subfuncs)
                if strcmp(env.subfuncs{k}.getSubFuncVal(),func)
                    paramSyms = env.subfuncs{k}.getSubFuncParamSyms();
                    paramVals = env.subfuncs{k}.getSubFuncParamVals();
                    for l=1:1:length(paramSyms)
                        func = regexprep(func,paramSyms(l),string(paramVals(l)));
                    end
                end
            end
        end
    end

    methods (Static)
        % env: Environment class ref
        function paramNames = getParamNames()
            paramNames = {'Model Runtime (t)','Incident Light (I)','Temperature (T)','Culture Volume (V)','Culture Surface Area (SA)';'t','I','T','V','SA'}';
        end
    end
end