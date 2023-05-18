classdef ODESys < handle
    properties
        species = struct(); % struct of Species in ODE system, (key: string, val: {Species})
        chemicals = struct(); % struct of Chemicals in ODE system, (key: string, val: {Chemical})

        environs = struct(); % struct of Environments in ODE system
        activeEnv = 'Photo_Bioreactor_PBR'; % string of field name of active Environment in environs struct
            % when pulling data for plotting, use this property to
            % reference the active environment with
            % sys.environs.(activeEnv)

        currentCompName = ""; % used to track changes in component/function editing
        currentFuncName = ""; % used to track changes in component/function editing
        currentFuncVal = ""; % used to track changes in component/function editing
        currentFuncCombo = ""; % used to track changes in component/function editing
        currentFuncType = ""; % used to track changes in component/function editing

        editing = false; % used for tracking changes in component/function editing

        degree = 0; % number of functions within the ODE system
        dydt = ""; % ODE system functions
        param = []; % ODE system parameter values
        f = {}; % ODE system subfuncs for passing into ode45
        helperFuncs = {}; % ODE system helper functions
        matches = []; % used to find parameters that are being regressed

        sysVars = {}; % system variable names
        plots = {}; % plot objects

        importedDataPath = ""; % path to user imported data
        importedData = {}; % user-imported data
        regParamList = {}; % list of parameters for easy access
        matchedVarsList = {}; % list of parameters that are matched in regression
        regSpecs = statset; % struct for regression specifications
        importedDataIdx = []; % for picking out which data to match in nlinfit
        regData = []; % storing regressed data
        reg_param_ct = 1;
        reg_analytics = struct( ...
                'beta',[], ...
                'R', [], ...
                'J', [], ...
                'CovB', [], ...
                'MSE', 0, ...
                'ErrorModelInfo', [] ...
            );
    end

    properties(Constant)
        sysFuncDef = containers.Map('KeyType','single','ValueType','char'); % containers.Map to track what Species or Chemical each function in dydt represents
    end

    methods
        function sys = ODESys()
            disp('ODESys init');

            % set default values
            sys.environs.Photo_Bioreactor_PBR = Environment("Photo_Bioreactor_PBR",sys.getModelVarNames());
            sys.environs.Eutrophic_Lake = Environment("Eutrophic_Lake",sys.getModelVarNames());
            sys.environs.Test_Tube = Environment("Test_Tube",sys.getModelVarNames());
            sys.environs.Shaking_Flask = Environment("Shaking_Flask",sys.getModelVarNames());
            sys.environs.Outdoor_Tank = Environment("Outdoor_Tank",sys.getModelVarNames());

            sys.regSpecs.paramIGs = [];
        end

        % sys: ODSSys class ref, key: string, val: cell array of
        %   conditions + units
        function updateInitCond(sys, key, val)
            sys.initCond(key) = val;
        end

        % sys: ODESys class ref, num: number, name: string, initConc:
        % number
        function sys = addSpecies(sys, num, name, initConc)
            % increases the degree of the ODE system
            sys.degree = sys.degree + 1;
            % maps Species to number in ODE system
            sys.sysFuncDef(sys.degree) = name;
            % instantiates Species, stores in ODESys Map
            specName = regexprep(name, ' ', '_');
            sys.species.(specName) = Component(name, num, initConc, "spec");
        end

        % sys: ODESys class ref, num: number, name: string
        function sys = addChemical(sys, num, name, initConc)
            % increases the degree of the ODE system
            sys.degree = sys.degree + 1;
            % maps Chemical to number in ODE system
            sys.sysFuncDef(sys.degree) = name;
            % instantiates Chemical, stores in ODESys Map
            chemName = replace(regexprep(regexprep(replace(regexprep(name,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            sys.chemicals.(chemName) = Component(name, num, initConc, "chem");
        end

        % sys: ODESys class ref, name: string
        function sys = addEnvironment(sys, name)
            % creates Environment field name, adds from ODESys Map
            envName = regexprep(regexprep(regexprep(regexprep(name, ' ', '_'), '(',''), ')',''), '-','_');
            % adds new Environment
            sys.environs.(envName) = Environment(name,sys.getModelVarNames());
        end

        % sys: ODESys class ref, name: string, field: string, val: number | string
        function sys = updateSpecies(sys, name, field, val)
            specName = regexprep(name, ' ', '_');
            spec = sys.species.(specName);
            switch field
                case "initConc"
                    spec.setInitConc(val);
                    
            end
        end

        % sys: ODESys class ref, name: string, field: string, val: number | string
        function sys = updateChemical(sys, name, field, val)
            chemName = replace(regexprep(regexprep(replace(regexprep(name,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            chem = sys.chemicals.(chemName);
            switch field
                case "initConc"
                    chem.setInitConc(val);
                case "nutModel"
                    chem.setNutModel(val);
            end
        end

        % sys: ODESys class ref, name: string, field: string, val: number |
        % string
        function updateEnvironment(sys, name, field, val)
            envName = regexprep(regexprep(regexprep(regexprep(name, ' ', '_'), '(',''), ')',''), '-','_');
            env = sys.environs.(envName);
            switch field
                case "Incident Light (I)"
                    env.setLightFunc(val);
                case "Temperature (T)"
                    env.setTempFunc(val);
                case "Initial Culture Volume (V)"
                    env.setCulVol(val);
                case "Culture Surface Area (SA)"
                    env.setCulSA(val);
                case "Model Time (t)"
                    env.setModelTime(val);
            end

            sys.updateEnvSubFuncs(val,field);
        end

        % sys: ODESys class ref, name: string
        function sys = removeSpecies(sys, name)
            % creates Species field name, removes from ODESys Map
            specName = regexprep(name, ' ', '_');
            sys.species = rmfield(sys.species, specName);
        end

        % sys: ODESys class ref, name: string
        function sys = removeChemical(sys, name)
            % creates Chemical field name, removes from ODESys Map
            chemName = replace(regexprep(regexprep(replace(regexprep(name,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            sys.chemicals = rmfield(sys.chemicals, chemName);
        end

        % sys: ODESys class ref, name: string
        function sys = removeEnvironment(sys, name)
            % creates Environment field name, removes from ODESys Map
            envName = regexprep(regexprep(regexprep(regexprep(name, ' ', '_'), '(',''), ')',''), '-','_');
            sys.environs = rmfield(sys.environs, envName);
        end

        % sys: ODESys class ref, name: string
        function sys = setEnvironment(sys, name)
            % creates Environment field name
            envName = regexprep(regexprep(regexprep(regexprep(name, ' ', '_'), '(',''), ')',''), '-','_');
            sys.activeEnv = envName;

            % ### FIXME: add code here to update the EnvSubFuncs based on
            % active environment
        end

        % sys: ODESys class ref
        function [paramVals, paramUnits, funcHandles] = getCurrentEnvironParams(sys)
            paramVals = sys.environs.(sys.activeEnv).getParamVals();
            funcHandles = sys.environs.(sys.activeEnv).getParamFuncHandles(sys.getModelVarNames());
            custDefault = EnvDefaults();
            units = custDefault.Units;
            paramUnits = {units.modelTime;units.lightFunc;units.tempFunc;units.culVol;units.culSA};
        end

        % sys: ODESys class ref, prop: string
        function res = getSpecies(sys, prop)
            switch prop
                case 'name'
                    specs = struct2cell(sys.species);
                    % need to write loop for this
                    res = cell(length(specs),1);
                    for i=1:1:length(specs)
                        res{i} = specs{i}.name;
                    end
                case 'comp'
                    res = struct2cell(sys.species);
            end
        end

        % sys: ODESys class ref, prop: string
        function res = getChemicals(sys, prop)
            switch prop
                case 'name'
                    chems = struct2cell(sys.chemicals);
                    % need to write loop for this
                    res = cell(length(chems),1);
                    for i=1:1:length(chems)
                        res{i} = chems{i}.name;
                    end
                case 'comp'
                    res = struct2cell(sys.chemicals);
            end
        end

        % sys: ODESys class ref, modelName: string, chem: string,
        % interact: string, spec: string, func: string, funcName: string
        function updateModel(sys, compName, funcVal, funcName, funcCombo, funcType)
            % get component object
            comp = sys.getCompByName(compName);

            % gets the param expressions of dependent variables
            depVars = sys.getModelVarNames();

            % param and depVar separation for funcVal
            symVarStrArr = string(symvar(str2sym(funcVal)));
            % creates array of system variables as strings
            varStr = string(intersect(symVarStrArr,depVars));
            % creates array of parameters as strings
            paramStr = string(setxor(symVarStrArr,varStr));

            % adding model to Comp
            comp.setModel(funcVal, funcName, funcCombo, funcType, paramStr, varStr, sys.editing, sys.currentFuncName);

            % interpret which terms are variables in model function
%             vars = convertCharsToStrings(split(func,["+","-","*","/","^","(",")"]));
%             varChars = {};
%             paramChars = {};
%             j = 1;
%             k = 1;
%             for i=1:1:length(vars)
%                 % do we need both 1st and 3rd logic gate here?
%                 if isnan(str2double(vars(i))) && vars(i) ~= ""
%                     if any(strcmp(convertCharsToStrings(depVars(:,2)),vars(i)))
%                         varChars{1,j} = vars(i); %#ok<AGROW> 
%                         j = j+1;
%                     else
%                         paramChars{1,k} = vars(i); %#ok<AGROW>
%                         k = k+1;
%                     end
%                 end
%             end
            % removing duplicates
%             varStr = unique([varChars{:}]);
%             paramStr = unique([paramChars{:}]); % use this for displaying param expressions later

            % converting function string to function_handle
%             funcStr = '@(';
%             for i=1:1:length(sys.tempParamStore)
%                 funcStr = [funcStr,app.tempParamStore{i}]; %#ok<AGROW> 
%             end
%             funcStr = [funcStr,')',func];

            % returning LaTeX output for UI - should be already done from
            % the SubFuncVF value changing callback
%             FFuncStr = "$F(";
%             GFuncStr = "$G(";
%             for i=1:1:length(FvarStr)
%                 if i==length(FvarStr)
%                     FFuncStr = FFuncStr+FvarStr(i);
%                 elseif i~=length(FvarStr)
%                     FFuncStr = FFuncStr+FvarStr(i)+",";
%                 end
%             end
%             for i=1:1:length(GvarStr)
%                 if i==length(GvarStr)
%                     GFuncStr = GFuncStr+GvarStr(i);
%                 elseif i~=length(GvarStr)
%                     GFuncStr = GFuncStr+GvarStr(i)+",";
%                 end
%             end
%             % converts function string into LaTeX
%             FFunc = sys.convertToLaTeX(FFunc);
%             GFunc = sys.convertToLaTeX(GFunc);
%             FFuncStr = FFuncStr+")="+FFunc+"$";
%             GFuncStr = GFuncStr+")="+GFunc+"$";
        end

        % sys: ODESys class ref, funcName: string
        function removeModel(sys)
            % get comp by name
            comp = sys.getCompByName(sys.currentCompName);

            % remove model
            comp.removeModel(sys.currentFuncVal,sys.currentFuncName,sys.currentFuncCombo,sys.currentFuncType);
        end

        % sys: ODESys class ref, funcName: string, funcSym: string,
        % funcVal: string, add: boolean
        function addRmHelperFuncs(sys,funcName,funcSym,funcVal,add)
            if add
                sys.helperFuncs{end+1} = SubFunc(funcVal,funcName,funcSym,0,0);
                sys.helperFuncs{end}.initParams(sys.getModelVarNames());
            else
                for k=1:1:length(sys.helperFuncs)
                    if strcmp(sys.helperFuncs{k}.getSubFuncName(),funcName) && ... 
                            strcmp(sys.helperFuncs{k}.getSubFuncSym(),funcSym)
                        sys.helperFuncs(k) = [];
                        break;
                    end
                end
            end
        end
        
        % sys: ODESys class ref, prop: string, funcName: string,
        % funcSym: string,  funcVal: string
        function updateHelperFuncs(sys,prop,funcName,funcSym,funcVal)
            switch prop
                case "name"
                    for k=1:1:length(sys.helperFuncs)
                        if strcmp(sys.helperFuncs{k}.getSubFuncSym(),funcSym) && ...
                                strcmp(sys.helperFuncs{k}.getSubFuncVal(),funcVal)
                            sys.helperFuncs{k}.setSubFuncName(funcName);
                        end
                    end
                case "sym"
                    for k=1:1:length(sys.helperFuncs)
                        if strcmp(sys.helperFuncs{k}.getSubFuncName(),funcName) && ...
                                strcmp(sys.helperFuncs{k}.getSubFuncVal(),funcVal)
                            sys.helperFuncs{k}.setSubFuncSym(funcSym);
                        end
                    end
                case "val"
                    for k=1:1:length(sys.helperFuncs)
                        if strcmp(sys.helperFuncs{k}.getSubFuncName(),funcName) && ...
                                strcmp(sys.helperFuncs{k}.getSubFuncSym(),funcSym)
                            sys.helperFuncs{k}.setSubFuncVal(funcVal);
                        end
                    end
            end
        end

        % sys: ODESys class ref, funcName: string
        function [funcSym,funcVal] = getHelperFunc(sys,funcName)
            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName,funcName)
                    funcSym = sys.helperFunc{k}.getSubFuncSym();
                    funcVal = sys.helperFunc{k}.getSubFuncVal();
                    break;
                end
            end
        end

        % sys: ODESys class ref, compName: string
        function funcs = getCompFuncNames(sys,compName,funcCombo)
            % get comp name
            comp = sys.getCompByName(compName);

            % get funcs from comp
            funcs = comp.getFuncNames(funcCombo);
        end

        % sys: ODESys class ref, compName: string
        function sym = getCompSymByName(sys,compName)
            comp = sys.getCompByName(compName);
            sym = comp.getSym();
        end

        % sys: ODESys class ref, compName: string, funcName: string
        function funcVal = getCompFuncByFuncName(sys,compName,funcName)
            comp = sys.getCompByName(compName);
            funcVal = comp.getFuncValByName(funcName);
        end

        % sys: ODEsys class ref, compName: string, valName: string,
        % editing: boolean, editedFunc: string
        function existFunc = getGovFuncByCompName(sys,compName,funcCombo,valsOrNames)
            comp = sys.getCompByName(compName);
            funcs = comp.getGovFuncs(funcCombo,valsOrNames,sys.editing,sys.currentFuncVal,sys.currentFuncName);
            existFunc = strings(1,length(funcs));
            for k=1:1:length(funcs)
                if funcCombo == "Multiplicative"
                    for l=1:1:length(funcs{k})
                        existFunc(k) = existFunc(k)+"\left("+latex(str2sym(funcs{k}{l}))+"\right)";
                    end
                elseif funcCombo == "Additive"
                    for l=1:1:length(funcs{k})
                        if l == 1
                            existFunc(k) = existFunc(k)+"\left("+latex(str2sym(funcs{k}{l}))+"\right)";
                        else
                            existFunc(k) = existFunc(k)+"+\left("+latex(str2sym(funcs{k}{l}))+"\right)";
                        end
                    end
                end
            end
        end

        % sys: ODESys class ref, compName: string, funcCombo: string,
        % funcType: string
        function funcName = getDefaultFuncName(sys,compName)
            comp = sys.getCompByName(compName);
            funcNum = comp.getFuncCount();
            funcName = "F"+funcNum;
        end

        % sys: ODESys class ref, compName: string, funcCombo: string,
        % funcType: string
        function funcVal = getDefaultFuncVal(sys,compName,funcCombo,funcType)
            comp = sys.getCompByName(compName);
            funcVal = CompDefaults.getDefaultFuncVals(comp,funcCombo,funcType);
        end

        % sys: ODESys class ref, compName: string, funcCombo: string
        function funcs = getValidFuncs(sys,compName,funcCombo)
            comp = sys.getCompByName(compName);
            funcs = CompDefaults.getValidFuncs(comp,funcCombo);
        end

        % sys: ODESys class ref, currentFuncName:
        % string, currentFuncVal: string, currentFuncCombo: string,
        % currentFuncType: string
        function setCurrent(sys,currentCompName,currentFuncName,currentFuncVal,currentFuncCombo,currentFuncType)
            if currentFuncName ~= ""
                sys.currentCompName = currentCompName;
                sys.currentFuncName = currentFuncName;
                sys.currentFuncVal = currentFuncVal;
                sys.currentFuncCombo = currentFuncCombo;
                sys.currentFuncType = currentFuncType;
            end
        end

        % sys: ODESys class ref
        function [currentCompName,currentFuncName,currentFuncVal,currentFuncCombo,currentFuncType] = getCurrent(sys)
            currentCompName = sys.currentCompName;
            currentFuncName = sys.currentFuncName;
            currentFuncVal = sys.currentFuncVal;
            currentFuncCombo = sys.currentFuncCombo;
            currentFuncType = sys.currentFuncType;
        end

        % sys: ODESys class ref, edited: boolean
        function setEditing(sys,editing)
            sys.editing = editing;
        end

        % sys: ODESys class ref
        function editing = getEditing(sys)
            editing = sys.editing;
        end

        % sys: ODESys class ref, paramSym: string,
        % newVal: number, newUnit: string, newParamName: string, funcName:
        % string
        function updateMultiParams(sys,paramSym,newVal,newUnit,newParamName,funcName)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            for k=1:1:length(comps)
                grthParams = comps{k}.getGrthParams();
                if any(strcmp(grthParams(:,3),paramSym))
                    comps{k}.updateParam(paramSym,newVal,newUnit,newParamName,funcName);
                end
            end

            envs = struct2cell(sys.environs);
            for k=1:1:length(envs)
                for l=1:1:length(envs{k}.subfuncs)
                    envFuncName = envs{k}.subfuncs{l}.getSubFuncName();
                    grthParams = envs{k}.getGrthParamsByEnvFuncName(envFuncName);
                    if any(strcmp(grthParams(:,3),paramSym))
                        envs{k}.updateEnvFuncParams(envFuncName,paramSym,newVal,newUnit,newParamName);
                    end
                end
            end

            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName(),helperFuncName)
                    helperFuncName = sys.helperFncs{k}.getSubFuncName();
                    grthParams = sys.getGrthParamsByHelperFuncName(helperFuncName);
                    if any(strcmp(grthParams(:,3),paramSym))
                        sys.helperFuncs{k}.updateParams(newParamName,paramSym,newVal,newUnit);
                    end
                end
            end
        end

        % sys: ODESys class ref, prop: char, val: string
        function type = getComponentType(sys, prop, val)
            % adjusting names for proper syntax
%             chemName = replace(regexprep(regexprep(replace(regexprep(val,' ','_'),'^','_'),'+','p'),'-','n'),'.','_')
%             specName = regexprep(val, ' ', '_')
            % gets all species and chemical names
            specs = sys.getSpecies(prop);
            chems = sys.getChemicals(prop);
            if any(strcmp(specs,val))
                type = "spec";
            elseif any(strcmp(chems,val))
                type = "chem";
            else
                type = "none";
            end
        end

        % sys: ODESys class ref, compName: string
        function comp = getCompByName(sys,compName)
            compType = sys.getComponentType('name',compName);
            if compType == "spec"
                specName = regexprep(compName, ' ', '_');
                comp = sys.species.(specName);
            elseif compType == "chem"
                chemName = replace(regexprep(regexprep(replace(regexprep(compName,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
                comp = sys.chemicals.(chemName);
            else
                comp = [];
                disp("Doesn't Work, ODESys.m: line 295");
                return;
            end
        end

        % sys: ODESys class ref, type: string
        function names = getCompNamesByType(sys,type)
            if type == "spec"
                names = sys.getSpecies('name');
            elseif type == "chem"
                names = sys.getChemicals('name');
            end
        end

        % sys: ODESys class ref
        function names = getEnvFuncNames(sys)
            names = sys.environs.(sys.activeEnv).getParamNames();
        end

        % sys: ODESys class ref
        function names = getHelperFuncNames(sys)
            names = string.empty(0,1);
            for k=1:1:length(sys.helperFuncs)
                names(k,1) = sys.helperFuncs{k}.getSubFuncName();
            end
        end

        % sys: ODESys class ref
        % Model Param Nomenclature:
        % Species concentrations: X1, X2, ... , Xn
        % Chemical concentrations: C1, C2, ... , Cn
        % Light: I
        % Time: t
        % Temperature: T
        % Environment Volume: V
        % Environment Surface Area: SA
        function vars = getModelVarNames(sys,varargin)
            % compiling all parameters from Environment and Species and
            % Chemical concentrations to push into Model Param Table
            specs = struct2cell(sys.species);
            chems = struct2cell(sys.chemicals);

            onlySpecsChem = false;
            if ~isempty(varargin)
                onlySpecsChem = varargin{1};
                vars = cell(length(specs)+length(chems),2); % 5 is from # of Environment properties
            else
                vars = cell(length(specs)+length(chems)+length(Environment.getParamNames()),2); % 5 is from # of Environment properties
            end

            % Species Concentrations
            if length(specs) > 0 %#ok<ISMT> 
                for i=1:1:length(specs)
                    vars{i,1} = specs{i}.name;
                    vars{i,2} = ['X',num2str(i)];
                end
            end
            % Chemical Concentrations
            if length(chems) > 0 %#ok<ISMT> 
                for i=1:1:length(chems)
                    vars{i+length(specs),1} = chems{i}.name;
                    vars{i+length(specs),2} = ['C',num2str(i)];
                end
            end
            if ~onlySpecsChem
                % Environmental parameters
                envParams = Environment.getParamNames();
                for i=1:1:size(envParams,1)
                    vars{i+length(specs)+length(chems),1} = envParams{i,1};
                    vars{i+length(specs)+length(chems),2} = envParams{i,2};
                end

                for k=1:1:length(sys.helperFuncs)
                    vars{k+length(specs)+length(chems)+size(envParams,1),1} = sys.helperFuncs{k}.getSubFuncName();
                    vars{k+length(specs)+length(chems)+size(envParams,1),2} = sys.helperFuncs{k}.getSubFuncSym();
                end
            end
        end

        % sys: ODESys class ref, name: string
        function params = getGrthParamsByCompName(sys,name)
            comp = sys.getCompByName(name);
            params = comp.getGrthParams();
        end

        % sys: ODESys class ref, envFuncName: string
        function params = getGrthParamsByEnvFuncName(sys,envFuncName)
            params = sys.environs.(sys.activeEnv).getGrthParamsByEnvFuncName(envFuncName);
        end

        % sys: ODESys class ref, helperFuncName: string
        function params = getGrthParamsByHelperFuncName(sys,helperFuncName)
            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName(),helperFuncName)
                    paramNames = sys.helperFuncs{k}.getSubFuncParamNames();
                    paramSyms = sys.helperFuncs{k}.getSubFuncParamSyms();
                    paramVals = sys.helperFuncs{k}.getSubFuncParamVals();
                    paramUnits = sys.helperFuncs{k}.getSubFuncParamUnits();
                    paramNums = 1:1:length(paramNames);
                    params = cell(length(paramNums),6);
                    for l=paramNums
                        params{l,1} = char(string(paramNums(l)));
                        params{l,2} = helperFuncName;
                        params{l,3} = char(paramSyms(l));
                        params{l,4} = char(paramNames(l));
                        params{l,5} = paramVals(l);
                        params{l,6} = char(paramUnits(l));
                    end
                end
            end
        end

        % sys: ODESys class ref
        function params = getAllGrthParams(sys)
            params = {};
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            for k=1:1:length(comps)
                for l=1:1:length(comps{k}.funcParams)
                    for m=1:1:length(comps{k}.funcParams{l}.params)
                        params{:,1} = [params{:,1},string(comps{k}.funcParams{l}.params{m}.sym)];
                        params{:,2} = [params{:,2},comps{k}.funcParams{l}.params{m}.val];
                    end
                end
            end
        end

        % sys: ODESys class ref, specName: string
        function specNum = getSpeciesNumFromName(sys,specName)
            name = regexprep(specName, ' ', '_');
            specNum = sys.species.(name).number;
        end

        % sys: ODESys class ref, chemName: string
        function chemNum = getChemicalNumFromName(sys, chemName)
            name = replace(regexprep(regexprep(replace(regexprep(chemName,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            chemNum = sys.chemicals.(name).number;
        end

        % sys: ODESys class ref, specNum: number
        function specName = getSpeciesNameFromNum(sys, specNum)
            specs = struct2cell(sys.species);
            specName = specs{specNum}.name;
        end

        % sys: ODESys class ref, chemNum: number
        function chemName = getChemicalNameFromNum(sys, chemNum)
            chems = struct2cell(sys.chemicals);
            chemName = chems{chemNum}.name;
        end

        % ### FIXME: plot specifications
        % need to add code to allow user to specify a t range and precision
        % of t

        % sys: ODEsys class ref, tRange: number[], tPtsNb: number
        function [tRes,yRes] = compileModel(sys)
            % 1. loop through funcParams in each Component and assemble
            % functions into full governing function
            %   replace each parameter name with p({gloNum})
            %   can get the parameters to be converted from string to sym
            %   with all parameters saved in the p array by:
            %       str2sym(func+"p({gloNum}"+func)
            %   also need to convert names of system variables
            % 2. convert each governing function to sym function
            % 3. create ODESys of sym functions

            % ### FIXME: need to test to make sure it works
            % ### FIXME: need to explicitly add multiplication symbol
            % between multiplicative terms

            sys.dydt = "@(t,y,p,f) [";
            sys.param = [];
            sys.f = {};
            comps = [sys.getSpecies('comp');sys.getChemicals('comp')];
            sysVar = sys.getModelVarNames();
            for k=1:1:length(comps)
                % add logic to compile component gov funcs in order that
                % allows for growth-associated product funcs to be a
                % function of biomass funcs, and substrate funcs to be a
                % function of biomass and product funcs
                % ### FIXME: requires special logic/restrictions in the
                % governing function portion
                [govFunc, params, ~] = comps{k}.compileGovFunc(length(sys.param),sys.reg_param_ct,{},false);
                for l=1:1:length(sys.helperFuncs)
                    govFunc = regexprep(govFunc,sys.helperFuncs{l}.getSubFuncSym(),"f{"+(l+length(sys.environs.(sys.activeEnv).subfuncs))+"}(y,t,p,f)");
                end
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f)");
                    end
                end
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end

                if k ~= length(comps)
                    sys.dydt = sys.dydt + govFunc + ";";
                else
                    sys.dydt = sys.dydt + govFunc + "]";
                end
                sys.param = [sys.param,params];
            end

            % creating functions for environmental conditions
            env_param_ct = 1;
            for k=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                % ### FIXME: add feature to define helpers with another
                % helper (this is very complicated, for now just don't let
                % helpers be defined by helpers)
                govFunc = string(sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncVal());
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f)");
                    end
                end
                syms = sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncParamSyms();
                govFuncLength = strlength(govFunc);
                govFunc = split(govFunc,"#");
                for l=1:1:length(syms)
                    for m=1:1:length(govFunc)
                        if strlength(govFunc(m)) > 1 || govFuncLength == 1
                            govFunc(m) = regexprep(govFunc(m),syms(l),"p("+(length(sys.param)+env_param_ct)+")");
                        end
                    end
                    env_param_ct = env_param_ct + 1;
                end
                govFunc = strrep(strrep(strjoin(govFunc),"#","")," ","");
                funcArgText = "@(y,t,p,f)";
                sys.f{end+1} = str2func(funcArgText+govFunc);
                sys.param = [sys.param,sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncParamVals()];
            end

            % creating functions for helper functions
            helper_param_ct = 1;
            for k=1:1:length(sys.helperFuncs)
                % ### FIXME: add feature to define helpers with another
                % helper (this is very complicated, for now just don't let
                % helpers be defined by helpers)
                govFunc = string(sys.helperFuncs{k}.getSubFuncVal());
                for l=1:1:length(sys.helperFuncs)
                    govFunc = regexprep(govFunc,sys.helperFuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f)");
                end
                syms = sys.helperFuncs{k}.getSubFuncParamSyms();
                govFuncLength = strlength(govFunc);
                govFunc = split(govFunc,"#");
                for l=1:1:length(syms)
                    for m=1:1:length(govFunc)
                        if strlength(govFunc(m)) > 1 || govFuncLength == 1
                            govFunc(m) = regexprep(govFunc(m),syms(l),"p("+(length(sys.param)+helper_param_ct)+")");
                        end
                    end
                    helper_param_ct = helper_param_ct + 1;
                end
                govFunc = strrep(strrep(strjoin(govFunc),"#","")," ","");
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f)");
                    end
                end
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end
                funcArgText = "@(y,t,p,f)";
                sys.f{end+1} = str2func(funcArgText+govFunc);
                sys.param = [sys.param,sys.helperFuncs{k}.getSubFuncParamVals()];
            end

            % converting to function_handle
            sys.dydt = str2func(sys.dydt');

            % running model for each plot
            % ### FIXME: need to be able to plot all system variables
            % ### FIXME: include all system variables in ODE system?
            % ### FIXME: test the plotting functionality
            
            % ### FIXME: add feature to allow user to specify time
            % precision of model
            tEnd = sys.environs.(sys.activeEnv).getModelTime();
            tSmooth = linspace(0,tEnd,100);
            for k=1:1:length(sys.plots)
                plot_obj = sys.plots{k};
                axes = plot_obj.axes;
                if plot_obj.getPlotProp("display") == true || plot_obj.getPlotProp("download") == true
                    for l=1:1:length(axes)-1
                        if axes{k}.varIsIC
                            % run computation for each plot
                            y0_span_arr = {0,0};
                            y0_var_num = {0,0};
                            for m=1:1:length(axes)
                                if ~axes{m}.isDV
                                    y0_span_arr{m} = linspace(axes{m}.loEvalLim,axes{m}.upEvalLim,axes{m}.nbEvalPts);
                                    y0_var_num{m} = find(strcmp(sysVar,axes{m}.varNames));
                                end
                            end
        
                            % ### FIXME: need to add the specifications to control
                            % whether this is for var IC or for var values
                            res = cell(length(y0_span_arr{1}),length(y0_span_arr{2}));
                            for m=1:1:length(y0_span_arr{1})
                                for n=1:1:length(y0_span_arr{2})
                                    y0 = [];
                                    comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
                                    for o=1:1:(length(fieldnames(sys.species))+length(fieldnames(sys.chemicals)))
                                        if o == y0_var_num{1}
                                            y0(o) = y0_span_arr{1}{m}; %#ok<AGROW> 
                                        elseif o == y0_var_num{2}
                                            y0(o) = y0_span_arr{2}{n}; %#ok<AGROW> 
                                        else
                                            y0(o) = comps{o}.getInitConc(); %#ok<AGROW>
                                        end
                                    end
                                    [tRes,yRes] = sys.runModel(tSmooth,y0);
                                    res{m,n} = [yRes,tRes];
                                end
                            end
                        else
                            y0 = [];
                            comps = [sys.getSpecies('comp');sys.getChemicals('comp')];
                            for m=1:1:length(comps)
                                y0(m) = comps{m}.getInitConc(); %#ok<AGROW>
                            end
                            [tRes,yRes] = sys.runModel(tSmooth,y0);
                            res{1,1} = [yRes,tRes];
                        end
                    end

                    % plot models on fig
                    fig = figure('Name',plot_obj.title);
                    hold on;
                    if length(axes) == 2
                        xVarIdx = strcmp(sysVar(:,2),axes{1}.varNames);
                        yVarIdx = zeros(size(axes{2}.varNames));
                        for l=1:1:length(axes{2}.varNames)
                            yVarIdx(l) = find(strcmp(sysVar(:,2),axes{2}.varNames{l}));
                        end
                        plot(res{1,1}(:,xVarIdx),res{1,1}(:,yVarIdx));
                    elseif length(axes) == 3
                        % interpolation query points
                        qXVals = linspace(y0_span_arr{1}(1),y0_span_arr{1}(end),axes{1}.nbEvalPts);
                        qYVals = linspace(y0_span_arr{2}(1),y0_span_arr{2}(end),axes{2}.nbEvalPts)';
                        [qGridXVals,qGridYVals] = meshgrid(qXVals,qYVals);
                        
                        % 2D interpolation
                        DVVals = zeros(size(res));
                        DVIdx = strcmp(sysVar(2),{axes{3}.varNames});
                        for l=1:1:size(res,1)
                            for m=1:1:size(res,2)
                                DVVals(l,m) = res{l,m}(axes{3}.evaltVal,DVIdx);
                            end
                        end
                        interpZGrid = interp2(y0_span_arr{1}, y0_span_arr{2}, DVVals, qXVals, qYVals, 'makima');
                    
                        % plotting surface
                        surf(qGridXVals,qGridYVals,interpZGrid);
                    end
                    title(plot_obj.title);
                    for l=1:1:length(plot_obj.axes)
                        if l == 1
                            xlabel(axes{l}.title);
                            if ~axes{l}.useDefR
                                xlim(plot_obj.getDispLims(axes{l}));
                            end
                        elseif l == 2
                            ylabel(axes{l}.title);
                            if ~axes{l}.useDefR
                                ylim(plot_obj.getDispLims(axes{l}));
                            end
                        elseif l == 3
                            zlabel(axes{l}.title);
                            if ~axes{l}.useDefR
                                zlim(plot_obj.getDispLims(axes{l}));
                            end
                        end
                    end
                    legend(axes{end}.varNames);
                    hold off;

                    if plot_obj.getPlotProp("display") == false
                        close;
                    end
                    if plot_obj.getPlotProp("download") == true
                        plot_obj.downloadPlot(fig);
                    end
                end
            end
        end

        % sys: ODESys class ref, t: number[]
        function [tRes,yRes] = runModel(sys,t,y0)
            % iterating over BatchFunction with ode45
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [tRes,yRes] = ode45(@(t,y) sys.dydt(t,y,sys.param,sys.f),t,y0,opts);
        end

        % sys: ODESys class ref,
        function regStats = compileRegression(sys)
            sys.dydt = "@(t,y,p,f,r) [";
            sys.param = [];
            sys.reg_param_ct = 1;
            sys.f = {};
            comps = [sys.getSpecies('comp');sys.getChemicals('comp')];
            sysVar = sys.getModelVarNames();
            for k=1:1:length(comps)
                % add logic to compile component gov funcs in order that
                % allows for growth-associated product funcs to be a
                % function of biomass funcs, and substrate funcs to be a
                % function of biomass and product funcs
                % ### FIXME: requires special logic/restrictions in the
                % governing function portion
%                 idx = zeros(size(sys.regParamList,1),1);
%                 for l=1:1:size(sys.regParamList,1)
%                     if strcmp(sys.regParamList{l,2},comps{k}.name)
%                         idx(l) = 1;
%                     else
%                         idx(l) = 0;
%                     end
%                 end
%                 compRegParamList = sys.regParamList(find(idx),[1,6]);
                [govFunc, params, sys.reg_param_ct,regPUpdate] = comps{k}.compileGovFunc(length(sys.param),sys.reg_param_ct,sys.regParamList(:,[1,6]),true);
                if regPUpdate{1}, sys.regParamList{regPUpdate{1},6} = regPUpdate{2}; end
                for l=1:1:length(sys.helperFuncs)
                    govFunc = regexprep(govFunc,sys.helperFuncs{l}.getSubFuncSym(),"f{"+(l+length(sys.environs.(sys.activeEnv).subfuncs))+"}(y,t,p,f,r)");
                end
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f,r)");
                    end
                end
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end

                if k ~= length(comps)
                    sys.dydt = sys.dydt + govFunc + ";";
                else
                    sys.dydt = sys.dydt + govFunc + "]";
                end
                sys.param = [sys.param,params];
            end

            % creating functions for environmental conditions
            env_param_ct = 1;
            for k=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                % ### FIXME: add feature to define helpers with another
                % helper (this is very complicated, for now just don't let
                % helpers be defined by helpers)
                govFunc = string(sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncVal());
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f,r)");
                    end
                end
                syms = sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncParamSyms();
                govFuncLength = strlength(govFunc);
                govFunc = split(govFunc,"#");
                for l=1:1:length(syms)
                    for m=1:1:length(govFunc)
                        if strlength(govFunc(m)) > 1 || govFuncLength == 1
                            govFunc(m) = regexprep(govFunc(m),syms(l),"p("+(length(sys.param)+env_param_ct)+")");
                        end
                    end
                    env_param_ct = env_param_ct + 1;
                end
                govFunc = strrep(strrep(strjoin(govFunc),"#","")," ","");
                funcArgText = "@(y,t,p,f,r)";
                sys.f{end+1} = str2func(funcArgText+govFunc);
                sys.param = [sys.param,sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncParamVals()];
            end

            % creating functions for helper functions
            helper_param_ct = 1;
            for k=1:1:length(sys.helperFuncs)
                % ### FIXME: add feature to define helpers with another
                % helper (this is very complicated, for now just don't let
                % helpers be defined by helpers)
                govFunc = string(sys.helperFuncs{k}.getSubFuncVal());
                for l=1:1:length(sys.helperFuncs)
                    govFunc = regexprep(govFunc,sys.helperFuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f,r)");
                end
                syms = sys.helperFuncs{k}.getSubFuncParamSyms();
                govFuncLength = strlength(govFunc);
                govFunc = split(govFunc,"#");
                for l=1:1:length(syms)
                    for m=1:1:length(govFunc)
                        if strlength(govFunc(m)) > 1 || govFuncLength == 1
                            govFunc(m) = regexprep(govFunc(m),syms(l),"p("+(length(sys.param)+helper_param_ct)+")");
                        end
                    end
                    helper_param_ct = helper_param_ct + 1;
                end
                govFunc = strrep(strrep(strjoin(govFunc),"#","")," ","");
                for l=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                    if ~strcmp(sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,sys.environs.(sys.activeEnv).subfuncs{l}.getSubFuncSym(),"f{"+l+"}(y,t,p,f,r)");
                    end
                end
                for l=1:1:length(comps)
                    if l <= length(sys.getSpecies('comp'))
                        govFunc = regexprep(govFunc,"X"+l,"y("+l+")");
                    else
                        govFunc = regexprep(govFunc,"C"+(l-length(sys.getSpecies('comp'))),"y("+l+")");
                    end
                end
                funcArgText = "@(y,t,p,f,r)";
                sys.f{end+1} = str2func(funcArgText+govFunc);
                sys.param = [sys.param,sys.helperFuncs{k}.getSubFuncParamVals()];
            end

            % converting to function_handle
            sys.dydt = str2func(sys.dydt');
            sys.dydt = char(sys.dydt);
            sys.dydt(end) = ']';
            sys.dydt = string(sys.dydt);
            sys.dydt

            % converting to function_handle
            sys.dydt = str2func(sys.dydt');

            IVs = [sys.importedData{:,1}];
            sys.importedDataIdx = [];
            for k=1:1:length(sys.matchedVarsList)
                matchIdx = strcmp(cellstr(sysVar(:,2)),sys.matchedVarsList{k}.sysVarName);
                varSyms = cellstr(sysVar(:,2));
                if any(matchIdx) && ~strcmp("t",varSyms{matchIdx})
                    sys.importedDataIdx(end+1) = find(matchIdx);
                end
            end
            DVs = [sys.importedData{:,2:end}];
            % ### FIXME: upgrade later to allow automatic testing of
            % various starting guesses for each parameter
            beta0 = sys.regSpecs.paramIGs;

            y0 = [];
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            for m=1:1:(length(fieldnames(sys.species))+length(fieldnames(sys.chemicals)))
                y0(m) = comps{m}.getInitConc(); %#ok<AGROW>
            end
            [sys.reg_analytics.beta,sys.reg_analytics.R,sys.reg_analytics.J, ...
                sys.reg_analytics.CovB,sys.reg_analytics.MSE,sys.reg_analytics.ErrorModelInfo] = ...
                nlinfit(IVs,DVs,@(reg_param,t) sys.nLinRegHandler(reg_param,t,y0),beta0,sys.regSpecs);

            [tRes,yRes] = sys.runRegModel(IVs,y0,sys.reg_analytics.beta);
            sys.regData = [tRes,yRes];
            regStats = struct('importedDataIdx',sys.importedDataIdx, ...
                'beta',sys.reg_analytics.beta,'R',sys.reg_analytics.R, ...
                'J',sys.reg_analytics.J,'CovB',sys.reg_analytics.CovB,'MSE',sys.reg_analytics.MSE, ...
                'ErrorModelInfo',sys.reg_analytics.ErrorModelInfo, ...
                'IVs',IVs,'DVs',DVs,'tRes',tRes,'yRes',yRes);
        end

        % sys: ODESys class ref, param: number[], t: number[]
        function yRes = nLinRegHandler(sys,reg_param,t,y0)
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [~,yRes_all] = ode45(@(t,y) sys.dydt(t,y,sys.param,sys.f,reg_param),t,y0,opts);
            yRes = [];
            for k=1:1:length(sys.importedDataIdx)
                yRes(:,k) = yRes_all(:,sys.importedDataIdx(k)); %#ok<AGROW> 
            end
        end

        % sys: ODESys class ref, t: number[]
        function [tRes,yRes] = runRegModel(sys,t,y0,reg_param)
            % iterating over BatchFunction with ode45
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [tRes,yRes] = ode45(@(t,y) sys.dydt(t,y,sys.param,sys.f,reg_param),t,y0,opts);
        end

        % sys: ODESys class ref
        function reg_analytics = getRegAnalysis(sys)
            reg_analytics = sys.reg_analytics;
        end

        % sys: ODESys class ref
        function dydt = getODE(sys)
            dydt = sys.dydt;
        end

        % sys: ODESys class ref, funcVal: string, funcName: string
        function updateEnvSubFuncs(sys,funcVal,funcName)
            for k=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                if strcmp(sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncName(),funcName)
                    sys.environs.(sys.activeEnv).subfuncs{k}.setSubFuncVal(funcVal);
                    sys.environs.(sys.activeEnv).subfuncs{k}.findParamsInFunc(sys.getModelVarNames());
                end
            end
        end

        % sys: ODESys class ref, funcName: string, paramName: string,
        % paramSym: string, paramVal: number, paramUnit: string
        function updateEnvSubFuncParam(sys,funcName,paramName,paramSym,paramVal,paramUnit)
            for k=1:1:length(sys.environs.(sys.activeEnv).subfuncs)
                if sys.environs.(sys.activeEnv).subfuncs{k}.getSubFuncName() == funcName
                    sys.environs.(sys.activeEnv).subfuncs{k}.updateParams(paramName,paramSym,paramVal,paramUnit);
                end
            end
        end

        % sys: ODESys class ref, plotName: string
        function plot = getPlotByName(sys,plotName)
            for k=1:1:length(sys.plots)
                if strcmp(sys.plots{k}.title,plotName)
                    plot = sys.plots{k};
                    return;
                end
            end
        end

        % sys: ODESys class ref
        function [plot,axes] = createNewPlot(sys)
            sys.plots{end+1} = Plot("Plot "+(length(sys.plots)+1),sys.getModelVarNames(),sys.getModelVarNames(true));
            plot = sys.plots{end}.getAllPlotProps();
            axes = sys.plots{end}.getAllAxProps();
        end

        % sys: ODESys class ref
        function [plot,axes] = removePlot(sys,plotName,lastItem)
            sys.plots
            if lastItem
                sys.plots(1) = [];
                plot = {};
                axes = {};
            else
                for k=1:1:length(sys.plots)
                    if sys.plots{k}.title == plotName
                        sys.plots(k) = [];
                        break;
                    end
                end
                plot = sys.plots{1}.getAllPlotProps();
                axes = sys.plots{1}.getAllAxProps();
            end
        end

        % sys: ODESys class ref, plotName: string, title: string, axesNb: number
        function axes = updatePlotNameAx(sys,plotName,title,axesNb)
            plot = sys.getPlotByName(plotName);
            plot.updatePlot("title",title);
            if axesNb < plot.getPlotProp("axesNb")
                plot.removeZAxis();
            elseif axesNb > plot.getPlotProp("axesNb")
                plot.addZAxis();
            end
            plot.updatePlot("axesNb",axesNb);
            axes = plot.getAllAxProps();
        end

        % sys: ODESys class ref, plotName: string, dirName: string
        function setPlotDir(sys,plotName,dirName)
            plot = sys.getPlotByName(plotName);
            plot.updatePlot("downloadDir",dirName);
        end

        % sys: ODESys class ref,
        function setAllPlotDispDownloadFalse(sys)
            for k=1:1:length(sys.plots)
                sys.plots{k}.updatePlot("display",false);
                sys.plots{k}.updatePlot("download",false);
            end
        end

        % sys: ODESys class ref, plotName: string, display: boolean,
        % download: boolean
        function updatePlotDispDownload(sys,plotName,prop,val)
            plot = sys.getPlotByName(plotName);
            plot.updatePlot(prop,val);
        end
        
        % sys: ODESys class ref, plotName: string, axesDir: string, title:
        % string, varName: string, loDispLim: string | number, upDispLim:
        % string | number, useDefR: boolean
%         function updateAxProps(sys,plotName,axesDir,title,varName,loDispLim,upDispLim,useDefR)
%             plot = sys.getPlotByName(plotName);
%             if axesDir == "X"
%                 axesDir = 1;
%             elseif axesDir == "Y"
%                 axesDir = 2;
%             elseif axesDir == "Z"
%                 axesDir = 3;
%             end
%             props = ["title","varName","loDispLim","upDispLim","useDefR"];
%             vals = {title,varName,loDispLim,upDispLim,useDefR};
%             for k=1:1:length(props)
%                 plot.updateAx(axesDir,props(k),vals{k});
%             end
%         end

        % sys: ODESys class ref, plotName: string, axesDir: string, title:
        % string
        function axes = updateAxTitle(sys,plotName,axesDir,title)
            plot = sys.getPlotByName(plotName);
            if axesDir == "X"
                axesDir = 1;
            elseif axesDir == "Y"
                axesDir = 2;
            elseif axesDir == "Z"
                axesDir = 3;
            end
            plot.updateAx(axesDir,"title",title);
            axes = plot.getAllAxProps();
        end
        
        % sys: ODESys class ref, plotName: string, axesDir: string,
        % varName: string, addVar: boolean
        function varText = updateAxVars(sys,plotName,axesDir,varName,addVar)
            plot = sys.getPlotByName(plotName);
            splitRes = split(varName,'_');
            if length(splitRes) > 1
                varIsIC = true;
            else
                varIsIC = false;
            end

            if plot.getPlotProp("axesNb") == 2
                if axesDir == "X"
                    axesDir = 1;
                    plot.updateAx(axesDir,"varNames",varName);
                    plot.updateAx(axesDir,"varIsIC",varIsIC);
                elseif axesDir == "Y"
                    axesDir = 2;
                    if addVar
                        if ~any(strcmp(plot.getAxProp(axesDir,"varNames"),varName))
                            varNames = [plot.getAxProp(axesDir,"varNames"),{varName}];
                            plot.updateAx(axesDir,"varNames",varNames);
                        end
                    else
                        varNames = plot.getAxProp(axesDir,"varNames");
                        varNames(strcmp(varNames,varName)) = [];
                        plot.updateAx(axesDir,"varNames",varNames);
                    end
                end
            else
                if axesDir == "X"
                    axesDir = 1;
                    plot.updateAx(axesDir,"varNames",varName);
                    plot.updateAx(axesDir,"varIsIC",varIsIC);
                elseif axesDir == "Y"
                    axesDir = 2;
                    plot.updateAx(axesDir,"varNames",varName);
                    plot.updateAx(axesDir,"varIsIC",varIsIC);
                elseif axesDir == "Z"
                    axesDir = 3;
                    if addVar
                        if ~any(strcmp(plot.getAxProp(axesDir,"varNames"),varName))
                            varNames = [plot.getAxProp(axesDir,"varNames"),{varName}];
                            plot.updateAx(axesDir,"varNames",varNames);
                        end
                    else
                        varNames = plot.getAxProp(axesDir,"varNames");
                        varNames(strcmp(varNames,varName)) = [];
                        plot.updateAx(axesDir,"varNames",varNames);
                    end
                end
            end
            varText = "Variable(s): $";
            varNameStrings = convertCharsToStrings(plot.getAxProp(axesDir,"varNames"));
            for k=1:1:length(varNameStrings)
                if k == length(varNameStrings)
                    varText = varText + varNameStrings(k) + "$";
                    return;
                end
                varText = varText + varNameStrings(k) + ", ";
            end
        end

        % sys: ODESys class ref, plotName: string, axesDir: string,
        % loDispLim: string | number, upDispLim: string | number, useDefR:
        % boolean
        function updateAxDispLims(sys,plotName,axesDir,loDispLim,upDispLim,useDefR)
            plot = sys.getPlotByName(plotName);
            if axesDir == "X"
                axesDir = 1;
            elseif axesDir == "Y"
                axesDir = 2;
            elseif axesDir == "Z"
                axesDir = 3;
            end
            props = ["loDispLim","upDispLim","useDefR"];
            vals = {loDispLim,upDispLim,useDefR};
            for k=1:1:length(props)
                plot.updateAx(axesDir,props(k),vals{k});
            end
        end

        % sys: ODEsys class ref, plotName: string, axisName: string
        function axisData = getAxisICEvalData(sys,plotName,axisName)
            plot = sys.getPlotByName(plotName);
            axisData = plot.getAxisICEvalData(axisName);
        end

        % sys: ODESys class ref, plotName: string, axisName: string,
        % evaltVal: string | number, loEvalLim: string | number, upEvalLim:
        % string | number, nbEvalPts: number
        function setAxisICEvalData(sys,plotName,axisName,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            plot = sys.getPlotByName(plotName);
            plot.setAxisICEvalData(axisName,evaltVal,loEvalLim,upEvalLim,nbEvalPts);
        end

        % sys: ODESys class ref, plotName: string, axisName: string,
        % varName: string
        function varData = getVarICEvalData(sys,plotName,axisName,varName)
            plot = sys.getPlotByName(plotName);
            varData = plot.getVarICEvalData(axisName,varName);
        end

        % sys: ODESys class ref, plotName: string, axisName: string,
        % varName: string, evaltVal: string | number, loEvalLim: string | number, 
        % upEvalLim: string | number, nbEvalPts: number
        function setVarICEvalData(sys,plotName,axisName,varName,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            plot = sys.getPlotByName(plotName);
            plot.setVarICEvalData(axisName,varName,evaltVal,loEvalLim,upEvalLim,nbEvalPts);
        end

        % sys: ODESys class ref, filePath: string
        function setImportedData(sys,filePath,importedData)
            sys.importedDataPath = filePath;
            sys.importedData = importedData;
        end

        % sys: ODESys class ref
        function [importedData,path] = getUserImportedData(sys)
            importedData = sys.importedData;
            path = sys.importedDataPath;
        end

        % sys: ODESys class ref, compName: string, paramSym: string, updateType: string
        function params = updateRegParamList(sys,compName,paramSym,updateType)
            comp = sys.getCompByName(compName);
            paramList = comp.getGrthParams();
            cancel = false;
            if updateType == "Add"
                for k=1:1:size(sys.regParamList)
                    if strcmp(sys.regParamList{k,1},paramSym)
                        cancel = true;
                    end
                end
                
                if ~cancel
                    for k=1:1:size(paramList,1)
                        if strcmp(paramList{k,3},paramSym)
                            regParam = paramList(k,:);
                            break;
                        end
                    end
                    sys.regParamList{end+1,1} = regParam{3};
                    sys.regParamList{end,2} = compName;
                    sys.regParamList{end,3} = regParam{2};
                    sys.regParamList{end,4} = regParam{5};
                    sys.regParamList{end,5} = '~';
                    sys.regParamList{end,6} = "";

                    sys.regSpecs.DerivStep(end+1) = eps^(1/3);
                    sys.regSpecs.paramIGs(end+1) = regParam{5};
                end
            elseif updateType == "Remove"
                for k=1:1:size(sys.regParamList,1)
                    if strcmp(sys.regParamList{k,1},paramSym)
                        sys.regParamList(k,:) = [];
                        sys.regSpecs.DerivStep(k) = [];
                        sys.regSpecs.paramIGs(k) = [];
                        break;
                    end
                end
            end
            params = sys.regParamList(:,1:5);
        end

        % sys: ODESys, compName: string
        function params = removeRegParamByComp(sys,compName)
            for k=1:1:size(sys.regParamList,1)
                if strcmp(sys.regParamList{k,2},compName)
                    sys.regParamList(k,:) = [];
                    break;
                end
            end
            params = sys.regParamList;
        end

        % sys: ODESys class ref
        function regParamList = getRegParamList(sys)
            regParamList = sys.regParamList;
        end

        % sys: ODESys class ref, sysVar: string, importVar: string, importVarNum: num, match:
        % boolean
        function matchedVarsList = updateVarMatch(sys,sysVar,importVar,importVarNum,match)
            sysVarNames = string(sys.getModelVarNames());
            sysVar = string(sysVar);
            importVar = string(importVar);
            sysVarNum = 1;
            pairNum = 1;
            cancel = false;
            if match % only works if both variables aren't already taken
                for k=1:1:length(sys.matchedVarsList)
                    if sysVar == sys.matchedVarsList{k}.sysVarName || importVar == sys.matchedVarsList{k}.importVarName
                        cancel = true;
                        break;
                    end
                end
                if ~cancel
                    for k=1:1:length(sysVarNames)
                        if sysVarNames{k} == sysVar
                            sysVarNum = k;
                            break;
                        end
                    end
                    newPair = struct('sysVarName',sysVar,'sysVarNum',sysVarNum,'importVarName',importVar,'importVarNum',importVarNum);
                    sys.matchedVarsList{end+1} = newPair;
                end
            else
                for k=1:1:length(sys.matchedVarsList)
                    if sys.matchedVarsList{k}.sysVarName == sysVar && ...
                            sys.matchedVarsList{k}.importVarName == importVar
                        pairNum = k;
                        break;
                    end
                end
                sys.matchedVarsList(pairNum) = [];
            end
            matchedVarsList = cell(length(sys.matchedVarsList),2);
            for k=1:1:size(sys.matchedVarsList,2)
                matchedVarsList{k,1} = char(sys.matchedVarsList{k}.sysVarName);
                matchedVarsList{k,2} = char(sys.matchedVarsList{k}.importVarName);
            end
        end

        % sys: ODESys class ref, DVIdx: number, varName: string
        function [plotRegData,plotImportData] = updateRegPlotDV(sys,DVIdx)
            plotRegData = [];
            plotRegData(:,1) = sys.regData(:,1);
            plotRegData(:,2) = sys.regData(:,DVIdx+1);
            plotImportData = [];
            plotImportData(:,1) = [sys.importedData{:,1}];
            plotImportData(:,2) = [sys.importedData{:,sys.importedDataIdx(DVIdx)+1}];
        end

        % sys: ODESys class ref
        function regSpecs = getRegSpecs(sys)
            regSpecs = sys.regSpecs;
        end

        % sys: ODEsys class ref, regSpecs: statset
        function setRegSpecs(sys,regSpecs)
            sys.regSpecs = regSpecs;
        end
    end

    methods (Static)
        % funcAsStr: Boolean
        function [paramVals, paramUnits] = getDefaultEnvironParams(funcAsStr)
            custDefault = EnvDefaults();
            vals = custDefault.Custom_Env;
            units = custDefault.Units;
            disp({units.lightFunc;units.tempFunc;units.culVol;units.culSA;units.modelTime})
            paramVals = {vals.lightFunc;vals.tempFunc;vals.culVol;vals.culSA;vals.modelTime};
            if funcAsStr
                for i=1:1:length(paramVals)
                    if isa(paramVals{i},"function_handle")
                        paramVals{i} = func2str(paramVals{i});
                    end
                end
            end
            paramUnits = {units.lightFunc;units.tempFunc;units.culVol;units.culSA;units.modelTime};
        end

        % funcStr: string
        function res = convertToLaTeX(funcStr)
            res = latex(str2sym(funcStr));
        end

        % modelName: string, chemNum: number, specNum: number
        function [F,G] = getModelFunc(modelName, chemNum, specNum)
            if modelName == "Custom Model"
                % set ModelF value
                F = "1";
                G = "0";
            elseif modelName == "Monod"
                F = "(C"+chemNum+")/(K"+chemNum+"_"+specNum+"+C"+chemNum+")";
                G = "0";
            elseif modelName == "Linear"
                F = "A"+chemNum+"_"+specNum+"*C"+chemNum;
                G = "0";
            end
        end

        % sys: ODESys class ref, comp: appComponent, enabled: boolean
        function toggleChildrenEnabled(sys, comp, enabled)
            if length(allchild(comp.Children)) > 0 %#ok<ISMT> 
                for i=1:1:length(allchild(comp.Children))
                    sys.toggleChildrenEnabled(comp.Children(i), enabled);
                end
            end

            if isprop(comp, "Enable")
                comp.Enable = enabled;
            end
        end
    end
end