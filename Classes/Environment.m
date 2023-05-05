classdef Environment < handle
    properties
        name = '';

        subfuncs = {};
    end

    methods
        function env = Environment(name)
            env.name = name;

            env.setDefault(name);
        end

        % env: Environment class ref, val: number
        function env = setModelTime(env, val)
            env.subfuncs{1}.setSubFuncVal(val);
        end

        % env: Environment class ref, funcStr: string
        function env = setLightFunc(env, funcStr)
            env.subfuncs{2}.setSubFuncVal(funcStr);
        end

        % env: Environment class ref, funcStr: string
        function env = setTempFunc(env, funcStr)
            env.subfuncs{3}.setSubFuncVal(funcStr);
        end

        % env: Environment class ref, val: number
        function env = setCulVol(env, val)
            env.subfuncs{4}.setSubFuncVal(val);
        end

        % env: Environment class ref, val: number
        function env = setCulSA(env, val)
            env.subfuncs{5}.setSubFuncVal(val);
        end

        % env: Environment class ref, name: string
        function env = setDefault(env, name)
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
                env.subfuncs{k} = SubFunc(envParamVals{k},envParamNames{k,1},0,0);
                env.subfuncs{k}.updateParams(envParamNames{k,2},envParamNames{k,2},1,"");
            end
        end

        % env: Environment class ref
        function paramVals = getParamVals(env)
            paramVals = {env.subfuncs{1}.getSubFuncVal();env.subfuncs{2}.getSubFuncVal();env.subfuncs{3}.getSubFuncVal();env.subfuncs{4}.getSubFuncVal();env.subfuncs{5}.getSubFuncVal()};
        end
    end

    methods (Static)
        % env: Environment class ref
        function paramNames = getParamNames()
            paramNames = {'Model Runtime (t)','Incident Light (I)','Temperature (T)','Culture Volume (V)','Culture Surface Area (SA)';'t','I','T','V','SA'}';
        end
    end
end