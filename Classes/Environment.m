classdef Environment < handle
    properties
        name = '';

        lightFunc = ""; % function for incident light over time, ppf
        tempFunc = ""; % function for temperature over time, C
        culVol = ""; % volume of growth culture, L
        culSA = ""; % surface area of growth culture, m^2
        modelTime = ""; % time duration of model, hr
    end

    methods
        function env = Environment(name)
            env.name = name;

            env.setDefault(name);
        end
        % env: Environment class ref, funcStr: string
        function env = setLightFunc(env, funcStr)
            env.lightFunc = funcStr;
        end

        % env: Environment class ref, funcStr: string
        function env = setTempFunc(env, funcStr)
            env.tempFunc = funcStr;
        end

        % env: Environment class ref, val: number
        function env = setCulVol(env, val)
            env.culVol = val;
        end

        % env: Environment class ref, val: number
        function env = setCulSA(env, val)
            env.culSA = val;
        end

        % env: Environment class ref, val: number
        function env = setModelTime(env, val)
            env.modelTime = val;
        end

        % env: Environment class ref, name: string
        function env = setDefault(env, name)
            envDef = EnvDefaults();
            if ~any(envDef.Names == name)
                name = "Custom_Env";
            end
            env.setLightFunc(envDef.(name).lightFunc);
            env.setTempFunc(envDef.(name).tempFunc);
            env.setCulVol(envDef.(name).culVol);
            env.setCulSA(envDef.(name).culSA);
            env.setModelTime(envDef.(name).modelTime);
        end

        % env: Environment class ref
        function paramVals = getParamVals(env)
            paramVals = {env.lightFunc;env.tempFunc;env.culVol;env.culSA;env.modelTime};
        end
    end

    methods (Static)
        % env: Environment class ref
        function paramNames = getParamNames()
            paramNames = {{'Model Runtime (t)';'Incident Light (I)';'Temperature (T)';'Culture Volume (V)';'Culture Surface Area (SA)'}, ...
                {'t';'I';'T';'V';'SA'}};
        end
    end
end