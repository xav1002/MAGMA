classdef Environment < handle
    properties
        name = '';

        lightFunc = @(t) 1; % function for incident light over time, ppf
        tempFunc = @(t) 1; % function for temperature over time, C
        culVol = 1; % volume of growth culture, L
        culSA = 1; % surface area of growth culture, m^2
        modelTime = 1; % time duration of model, hr
    end

    methods
        function env = Environment(name)
            env.name = name;

            env.setDefault(name);
        end
        % env: Environment class ref, funcStr: string
        function env = setLightFunc(env, funcStr)
            % converts user's input string (without anonymous notation) to
            % function_handle
            if isa(funcStr,"string")
                funcStr = string(funcStr);
                func = str2func(strcat('@(t) ',funcStr));
            else
                func = funcStr;
            end
            env.lightFunc = func;
        end

        % env: Environment class ref, funcStr: string
        function env = setTempFunc(env, funcStr)
            % converts user's input string (without anonymous notation) to
            % function_handle
            if isa(funcStr,"string")
                funcStr = string(funcStr);
                func = str2func(strcat('@(t) ',funcStr));
            else
                func = funcStr;
            end
            env.tempFunc = func;
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
            disp(name)
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
            paramNames = {{'Time';'Incident Light';'Temperature';'Culture Volume';'Culture Surface Area'}, ...
                {'t';'I';'T';'V';'SA'}};
        end
    end
end