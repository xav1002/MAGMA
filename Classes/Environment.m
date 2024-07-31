classdef Environment < handle
    properties
        name = '';

        reactorType = "";
        reactorSpecificParams = {};

        subfuncs = {};

        T_comp = {};
        P_comp = {};
        V_comp = {};
        SV_comps = {};
        H3O_comp = {};
        OH_comp = {};

        K_w = 1E-14;
    end

    methods
        function env = Environment(name,modelVarSyms,defaultParamVals)
            env.name = name;

            env.setDefault(modelVarSyms,defaultParamVals);
        end

        % env: Environment class ref, funcStr: string
        function setLightFunc(env, funcStr)
            env.subfuncs{1}.setSubFuncVal(funcStr);
        end

        % env: Environment class ref, funcStr: string
        function setMaxVol(env,funcStr)
            env.subfuncs{2}.setSubFuncVal(char(string(funcStr)));
        end

        % env Environment class ref
        function comp = getTComp(env)
            comp = env.T_comp;
        end

        % env Environment class ref
        function comp = getPComp(env)
            comp = env.P_comp;
        end

        % env Environment class ref
        function comp = getVComp(env)
            comp = env.V_comp;
        end

        % env: Environment class ref
        function maxV = getMaxV(env)
            maxV = env.subfuncs{2};
        end

        % env: Environment class ref, prop: string
        function comp = getSVComps(env,prop)
            if strcmp(prop,'comp')
                comp = env.SV_comps;
            elseif strcmp(prop,'name')
                comp = [];
                for k=1:1:length(env.SV_comps)
                    comp(k) = env.SV_comps{k}.getName(); %#ok<AGROW>
                end
            end
        end

        % env: Environment class ref
        function comp = getH3OComp(env)
            comp = env.H3O_comp;
        end

        % env: Environment class ref
        function comp = getOHComp(env)
            comp = env.OH_comp;
        end

        % env: Environment class ref
        function comps = getAllEnvComps(env)
            comps = {env.getTComp(),env.getPComp(),env.getVComp(),env.getH3OComp(),env.getOHComp()};
        end

        % env: Environment class ref, initT: number
        function setInitT(env,initT)
            env.T_comp.setInitConc(initT,'K');
        end
        
        % env: Environment class ref, initP: number
        function setInitP(env,initP)
            env.P_comp.setInitConc(initP,'kPa');
        end

        % env: Environment class ref, initV: number
        function setInitV(env, initV)
            env.V_comp.setInitConc(initV,'L');
        end

        % env: Environment class ref, initpH: number, defaultParamVals: {}
        function setInitpH(env,initpH,defaultParamVals)
            H3O_initConc = 10.^(-initpH).*defaultParamVals.MW_H3O;
            OH_initConc = 1E-14./(H3O_initConc).*(defaultParamVals.MW_H3O./defaultParamVals.MW_OH);
            env.H3O_comp.setInitConc(H3O_initConc,'-');
            env.OH_comp.setInitConc(OH_initConc,'-');
        end

        % env: Environment class ref, solventName: string, initVol: num
        function addSVComp(env,solventName,initVol)
            % what's the sym + num on SV_comps?
            for k=1:1:length(env.SV_comps)
                if env.SV_comps{k}.getNum() ~= k
                    nextNum = k;
                end
            end
            newComp = Component(char(solventName),nextNum,initVol,'svol',false,modelVarSyms);
            env.SV_comps(nextNum+1:end+1) = env.SV_comps(nextNum:end);
            env.SV_comps{nextNum} = newComp;
        end

        % env: Environment class ref, solventName: string, initVol: num
        function updateSVComp(env,solventName,initVol)
            % update initial volume
            for k=1:1:length(env.SV_comps)
                if strcmp(env.SV_comps{k}.getName(),[char(solventName)])
                    env.SV_comps{k}.setInitConc(initVol,'L');
                    return;
                end
            end
        end

        % env: Environment class ref, solventName: string
        function removeSVComp(env,solventName)
            for k=1:1:length(env.SV_comps)
                if strcmp(env.SV_comps{k}.getName(),char(solventName))
                    env.SV_comps(k) = [];
                    return;
                end
            end
        end

        % env: Environment class ref
        function initCond = getInitCond(env)
            % Incident light
            % ### UPDATE: need to get initial condition for incident light
            initCond = {env.T_comp.getInitConc();env.P_comp.getInitConc();env.V_comp.getInitConc();-log10(double(string(env.H3O_comp.getInitConc())))};
        end

        % env: Environment class ref
        function setReactorSpecificParams(env,params)
            env.reactorType = params{1};
            env.reactorSpecificParams = params{2};
    
            % setting maximum volume value;
            env.subfuncs{2}.setSubFuncVal(char(string(env.reactorSpecificParams.maxVol)));
        end

        % env: Environment class ref
        function [reactorType,params] = getReactorSpecificParams(env)
            reactorType = env.reactorType;
            params = env.reactorSpecificParams;
        end

        % env: Environment class ref, envFuncName: string, paramSym: string,
        % newVal: number, newUnit: string, newParamName: string
        function updateEnvFuncParams(env,envFuncName,paramSym,newVal,newUnit,newParamName,defaultParamVals)
            for k=1:1:length(env.subfuncs)
                if strcmp(env.subfuncs{k}.getSubFuncName(),envFuncName)
                    env.subfuncs{k}.updateParams(newParamName,paramSym,newVal,newUnit,defaultParamVals);
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
                    paramEditable = env.subfuncs{k}.getSubFuncParamEditable();
                    paramNums = 1:1:length(paramNames);
                    params = cell(length(paramNums),6);
                    for l=paramNums
                        if strcmp(paramUnits(l),"mg/L")
                            val = paramVals(l).*1000;
                        elseif strcmp(paramUnits(l),"mg/g")
                            val = paramVals(l).*1000;
                        elseif strcmp(paramUnits(l),"g/mg")
                            val = paramVals(l).*0.001;
                        elseif strcmp(paramUnits(l),"cells/mg")
                            val = paramVals(l).*0.001;
                        elseif strcmp(paramUnits(l),"mg/cells")
                            val = paramVals(l).*1000;
                        else
                            val = paramVals(l);
                        end
                        params{l,1} = char(string(paramNums(l)));
                        params{l,2} = char(envFuncName);
                        params{l,3} = char(paramSyms(l));
                        params{l,4} = char(paramNames(l));
                        params{l,5} = val;
                        params{l,6} = char(paramUnits(l));
                        params{l,7} = char(string((paramEditable{l}==1)));
                    end
                end
            end
        end

        % env: Environment class ref, name: string, modelVarSyms: string{},
        % defaultParamVals: struct
        function setDefault(env, modelVarSyms, defaultParamVals)
            envDef = EnvDefaults();
            envParamVals = {};
            envParamVals{1} = envDef.Values.lightFunc;
            envParamVals{2} = envDef.Values.maxVol;
            envParamNames = env.getEnvSubfNames();
            for k=1:1:length(envParamVals)
                env.subfuncs{k} = SubFunc(envParamVals{k},envParamNames{1,k},envParamNames{2,k},0,0);
                env.subfuncs{k}.initParams(modelVarSyms,defaultParamVals);
            end

            env.reactorSpecificParams = {envDef.Values.PBRHt,envDef.Values.PBRWidth,envDef.Values.PBRDepth,envDef.Values.maxVol, ...
                envDef.Values.lgtAttnModel,envDef.Values.lgtBioAttn,envDef.Values.lgtChemAttn};

            % creating components for T, P, V
            env.T_comp = Component('Temperature',0,30,'C','temp',false,0,0,'',0,'',defaultParamVals);
            env.P_comp = Component('Pressure',0,1.01325,'kPa','press',false,0,0,'',0,'',defaultParamVals);
            env.P_comp.addModel("P/(V_m-V_tot)*(L_i-L_o)",'Pressure Dilution',[],'Main',defaultParamVals);
            env.V_comp = Component('Volume',0,char(string((100.*100.*5)./(1E3).*0.99)),'L','vol',false,0,0,'',0,'',defaultParamVals);

            % creating pH components
            env.H3O_comp = Component('Hydronium',0,1E-7,'M','acid',false,19.023,0,'',0,'',defaultParamVals);
            env.H3O_comp.addModel("k_pH*(H3O_eq-H3O)",'pH Equilibrium',"k_pH",'Main',defaultParamVals);
            env.OH_comp = Component('Hydroxide',0,1E-7,'M','base',false,17.007,0,'',0,'',defaultParamVals);
            env.OH_comp.addModel("k_pH*(OH_eq-OH)",'pH Equilibrium',"k_pH",'Main',defaultParamVals);
        end

        % env: Environment class ref
        function paramVals = getParamVals(env,defaultParamVals)
            paramVals = {};
            for k=1:1:length(env.subfuncs)-1
                func = env.subfuncs{k}.getSubFuncVal();
                paramVals{end+1,1} = char(string(func)); %#ok<AGROW>
            end
            env_comps = env.getAllEnvComps();
            for k=1:1:length(env_comps)-1
                if strcmp(env_comps{k}.getType(),'acid')
                    paramVals{end+1,1} = char(string(-log10(env_comps{k}.getInitConc()./defaultParamVals.MW_H3O))); %#ok<AGROW>
                else
                    paramVals{end+1,1} = char(string(env_comps{k}.getInitConc())); %#ok<AGROW>
                end
            end
            for k=1:1:length(env.reactorSpecificParams)
                paramVals{end+1,1} = env.reactorSpecificParams{k}; %#ok<AGROW>
            end
        end

        % env: Environment Class ref, varSyms: {}
        function funcHandles = getParamFuncHandles(env,varSyms,defaultParamVals)
            funcHandles = cell(size(env.subfuncs'));
            for k=1:1:length(funcHandles)
                func = env.subfuncs{k}.getSubFuncVal();
                env.subfuncs{k}.findParamsInFunc(varSyms,defaultParamVals);
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
            paramNames = {'Temperature','Pressure','Volume','Hydronium Concentration', ...
                'Hydroxide Concentration','Incident Light','Maximum Volume'; ...
                'T','P','V','H3O','OH','I_0','V_m'}';
        end

        function paramNames = getEnvSubfNames()
            paramNames = {'Incident Light','Maximum Volume','Average Light Path Length';'I_0','V_m','l_avg'};
        end
    end
end