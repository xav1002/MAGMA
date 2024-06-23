classdef ODESys < handle
    properties
        updateConsoleFunc = 1;

        save_type = "file";
        save_file_name = "";
        save_db_handle = "";
        save_model_name = "New Model";
        save_model_desc = "";
        save_model_bio_names = {};
        save_model_chem_names = {};
        save_model_ref_str = "";

        is_link = false; % used to determine whether this is an instance of MAGMA or MAGMALink UI

        species = struct(); % struct of Species in ODE system, (key: string, val: {Species})
        chemicals = struct(); % struct of Chemicals in ODE system, (key: string, val: {Chemical})

        environs = struct(); % struct of Environments in ODE system
        activeEnv = 'Photo_Bioreactor_PBR'; % string of field name of active Environment in environs struct
            % when pulling data for plotting, use this property to
            % reference the active environment with
            % sys.environs.(sys.activeEnv)
        model_runtime = 20; % hr

        currentCompName = ""; % used to track changes in component/function editing
        currentFuncName = ""; % used to track changes in component/function editing
        currentFuncVal = ""; % used to track changes in component/function editing
        currentFuncCombo = ""; % used to track changes in component/function editing
        currentFuncType = ""; % used to track changes in component/function editing

        editing = false; % used for tracking changes in component/function editing

        solver = ode;

        degree = 0; % number of functions within the ODE system
        dydt = ""; % ODE system functions
        param = []; % ODE system parameter values
        f = {}; % ODE system subfuncs for passing into ode15s
        helperFuncs = {}; % ODE system helper functions
%         matches = []; % used to find parameters that are being regressed
        y0 = [];
        v = []; % I/O values
        nonNegSysVars = {};
        solverNonNeg = [];

%         sysVars = {}; % system variable names
        plots = {}; % plot objects
        subplot_ct = 0;
        subplot_row_ct = [1];
        subplot_col_ct = [1];
        subplots = {};
        TLs = {};
        data_export_dir = "";
        data_export_tables = {};

        PBR_max_vol = 1; % L
        max_liq_vol = 1;
        environ_vars = {};
        var_in_key = {};
        var_out_key = {};
        S_nums = {};

        input_streams = {}; % used to keep track of streams being inputted into bioreactor, cell of InputStream objects
        output_streams = {}; % used to keep track of streams being outputed from bioreactor, cell of OutputStream objects
        has_liq_in = false;
        has_liq_out = false;
        has_gas_in = false;
        has_gas_out = false;
        inportCt = 0;
        outportCt = 0;
        IO_flowrates = struct();

        importedDataPath = ""; % path to user imported data
        importedData = {}; % user-imported data
        importedDataColNames = {}; % column names for user imported data

        regParamList = {}; % list of parameters for easy access
        matchedVarsList = {}; % list of parameters that are matched in regression
        regSpecs = struct;
        importedDataIdx = []; % for picking out which data to match in nlinfit
        regData = []; % storing regressed data
        IVs = [];
        DVs = [];
        mat_IVs = [];
        mat_DVs = [];
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

    methods
        function sys = ODESys(updateConsoleFunc)
            disp('ODESys init');

            sys.updateConsoleFunc = updateConsoleFunc;

            % set default values
            sys.environs.Photo_Bioreactor_PBR = Environment("Photo_Bioreactor_PBR",sys.getModelVarNames(),sys.getDefaultParamVals(true));
            sys.checkLiqMaxVol();

            % adding helper functions for pH equilibrium values and pH, pOH
            sys.addRmHelperFuncs('Equilibrium Hydrionium Concentration','H3O_eq',"(K_w/OH*(MW_OH)*(MW_H3O))",true);
            sys.addRmHelperFuncs('Equilibrium Hydroxide Concentration','OH_eq',"(K_w/H3O*(MW_H3O)*(MW_OH))",true);
            sys.addRmHelperFuncs('pH','pH',"-log10(H3O/MW_H3O)",true);
            sys.addRmHelperFuncs('pOH','pOH',"-log10(OH/MW_H3O)",true);
            sys.addRmHelperFuncs('Light Intensity','I',"I_0",true);

            % create V_tot helper
            sys.addRmHelperFuncs('Total Solvent Volume','V_tot',"V",true);

            sys.regSpecs.paramIGs = [];
            sys.regSpecs.params = {};
            sys.regSpecs.solver = 'nlinfit';
            sys.regSpecs.nlinfitSpecs = struct( ...
                'RobustWgtFun','', ...
                'Tune',[], ...
                'TolFun',1E-8, ...
                'TolX',1E-8, ...
                'MaxIter',100, ...
                'FunValCheck','off', ...
                'DerivStep',[] ...
            );
            sys.regSpecs.lsqcurvefitSpecs = struct( ...
                'paramBnds',[], ...
                'Algorithm','trust-region-reflective', ...
                'PlotFcns',[], ...
                'TolFun',1E-6, ...
                'TolX',1E-6, ...
                'MaxIter',1000, ...
                'MaxFunEvals',1000, ...
                'FunValCheck','off' ...
            );
            sys.regSpecs.fminsearchSpecs = struct( ...
                'PlotFcns',[], ...
                'TolFun',1E-4, ...
                'TolX',1E-4, ...
                'MaxIter',1000, ...
                'MaxFunEvals',1000, ...
                'FunValCheck','off' ...
            );
            sys.regSpecs.fminconSpecs = struct( ...
                'paramBnds',[], ...
                'Algorithm','interior-point', ...
                'PlotFcns',[], ...
                'TolFun',1E-6, ...
                'TolX',1E-6, ...
                'MaxIter',1000, ...
                'MaxFunEvals',1000, ...
                'FunValCheck','off' ...
            );
        end

        % sys: ODESys class ref, num: number, name: string, initConc:
        % number
        function sys = addSpecies(sys, num, name, initConc,initConcUnit)
            % increases the degree of the ODE system
            sys.degree = sys.degree + 1;
            % instantiates Species, stores in ODESys Map
            specName = regexprep(name, ' ', '_');
            % ### DEBUG: add logic here to make sure that "~" doesn't
            % appear in any names - throw error if user tries to use "~"
            % in names
            sys.species.(specName) = Component(name,num,initConc,initConcUnit,"Biological Solute",false,0,0,'',0,'',sys.getDefaultParamVals());
        end

        % sys: ODESys class ref, num: number, name: string
        function sys = addChemical(sys, num, name, initConc, initConcUnit, type, is_vol, MW, h_const, h_const_u, dh_const, dh_const_u)
            % increases the degree of the ODE system
            sys.degree = sys.degree + 1;
            % instantiates Chemical, stores in ODESys Map
            chemName = replace(regexprep(regexprep(replace(regexprep(name,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            
            sys.chemicals.(chemName) = Component(name, num, initConc, initConcUnit, type, is_vol, MW, h_const, h_const_u, dh_const, dh_const_u, sys.getDefaultParamVals());

            if strcmp(type,'Suspended Solid Sorbent')
                helperVal = sys.getHelperFunc('Total Solvent Volume');
                helperVal = helperVal + "+V_S_" + num;
                sys.updateHelperFuncs("val",'Total Solvent Volume','V_tot',helperVal);

                sys.updatePBRMaxVol();
            end
            
            % add a helper function for the h_const as a
            % function of temperature
            % need to have H_0 and A be parameters that cannot be
            % edited
            if is_vol
                HConstHelpers = sys.chemicals.(chemName).getHConstHelperFuncs();
                for k=1:1:length(HConstHelpers)
                    paramNames = HConstHelpers{k}.getSubFuncParamNames();
                    paramSyms = HConstHelpers{k}.getSubFuncParamSyms();
                    paramVals = HConstHelpers{k}.getSubFuncParamVals();
                    paramUnits = HConstHelpers{k}.getSubFuncParamUnits();
                    sys.addRmHelperFuncs(HConstHelpers{k}.funcName,HConstHelpers{k}.funcSym,HConstHelpers{k}.funcVal,true);
                    [~,~,helperIdx] = sys.getHelperFunc(HConstHelpers{k}.getSubFuncName());
                    for l=1:1:length(paramNames)
                        sys.helperFuncs{helperIdx}.updateParams(paramNames(l),paramSyms(l),paramVals(l),paramUnits(l),true,sys.getDefaultParamVals());
                    end
                end

                P_comp = sys.environs.(sys.activeEnv).getPComp();
                sys.addModel(P_comp.getName(), ...
                    "0",char(sys.chemicals.(chemName).sym+" liq_to_bulk"),'Main');

                % add and remove funcParam for liquid gas mass transfer
                sys.addModel([sys.chemicals.(chemName).getName(),' (Liquid Phase)'],"0",'Gas MT','Liquid');
                sys.addModel([sys.chemicals.(chemName).getName(),' (Gas Phase)'],"0",'Liquid MT','Gas');
            end

            % updates IOFuncs
            sys.updateIOFunc();
        end

        % sys: ODESys class ref, name: string, field: string, val: number | string
        function sys = updateSpecies(sys,compName,initConc,initConcUnit)
            compName = regexprep(compName, ' ', '_');
            comp = sys.species.(compName);
            comp.setInitConc(initConc,initConcUnit);
        end

        % sys: ODESys class ref, chemName: string, chemInitConc: num,
        % is_vol: boolean, h_const: num, h_const_u: string, dh_const: num,
        % dh_const_u: string
        function sys = updateChemical(sys,name,chemInitConc,initConcUnit,is_vol,MW,h_const,h_const_u,dh_const,dh_const_u)
            chemName = replace(regexprep(regexprep(replace(regexprep(name,' ','_'),'^','_'),'+','p'),'-','n'),'.','_');
            chem = sys.chemicals.(chemName);

            prev_is_vol = chem.is_vol;

            chem.updateComp(chemInitConc,initConcUnit,is_vol,MW,h_const,h_const_u,dh_const,dh_const_u,sys.getDefaultParamVals());

            sys.updatePBRMaxVol();

            if is_vol
                % add HConst helpers
                HConstHelpers = sys.chemicals.(chemName).getHConstHelperFuncs();
                for k=1:1:length(HConstHelpers)
                    paramNames = HConstHelpers{k}.getSubFuncParamNames();
                    paramSyms = HConstHelpers{k}.getSubFuncParamSyms();
                    paramVals = HConstHelpers{k}.getSubFuncParamVals();
                    paramUnits = HConstHelpers{k}.getSubFuncParamUnits();
                    if ~prev_is_vol, sys.addRmHelperFuncs(HConstHelpers{k}.funcName,HConstHelpers{k}.funcSym,HConstHelpers{k}.funcVal,true); end
                    [~,~,helperIdx] = sys.getHelperFunc(HConstHelpers{k}.getSubFuncName());
                    for l=1:1:length(paramNames)
                        sys.helperFuncs{helperIdx}.updateParams(paramNames(l),paramSyms(l),paramVals(l),paramUnits(l),true,sys.getDefaultParamVals());
                    end
                end

                if ~prev_is_vol
                    P_comp = sys.environs.(sys.activeEnv).getPComp();
                    sys.addModel(P_comp.getName(), ...
                        "0",char(sys.chemicals.(chemName).sym+" liq_to_bulk"),'Main');

                    % add and remove funcParam for liquid gas mass transfer
                    sys.addModel([sys.chemicals.(chemName).getName(),' (Liquid Phase)'],"0",'Gas MT','Liquid');
                    sys.addModel([sys.chemicals.(chemName).getName(),' (Gas Phase)'],"0",'Liquid MT','Gas');
                end
            elseif ~is_vol && prev_is_vol
                % remove HConst helpers
                HConstHelpers = chem.getHConstHelperFuncs();
                for k=1:1:length(HConstHelpers)
                    sys.addRmHelperFuncs(HConstHelpers{k}.funcName,HConstHelpers{k}.funcSym,HConstHelpers{k}.funcVal,false);
                end

                P_comp = sys.environs.(sys.activeEnv).getPComp();
                P_comp.removeModel("k_"+chem.sym+"_g*"+chem.h_const_sym+"*("+chem.sym+"-"+chem.bulk_gas_sym+")", ...
                    chem.sym+" liq_to_bulk",'Main');

                chem.removeModelByName('Gas MT','Liquid');
                chem.removeModelByName('Liquid MT','Gas');
            end

            % updates IOFuncs
            sys.updateIOFunc();
        end

        % sys: ODESys class ref, name: string, field: string, val: number |
        % string
        function newLgtMdlLaTeX = updateEnvironment(sys, field, val)
            env = sys.environs.(sys.activeEnv);
            newLgtMdlLaTeX = "$I=I_0$";
            switch field
                case "Incident Light (I_0(t))"
                    env.setLightFunc(val);
                case "Initial Temperature (T_0)"
                    env.setInitT(val);
                case "Initial Culture Volume (V_0)"
                    env.setInitV(val);
                    sys.updatePBRMaxVol();
                case "Initial Pressure (P_0)"
                    env.setInitP(val);
                case "Initial pH (pH_0)"
                    env.setInitpH(val,sys.getDefaultParamVals());
                case "Model Runtime (t)"
                    sys.model_runtime = double(string(val));
                case "reactorSpecificParams"
                    env.setReactorSpecificParams(val);
                    sys.updatePBRMaxVol();
                    newLgtMdlLaTeX = sys.updateLightAttenuationModel();
            end

            sys.updateEnvSubFuncs(val,field);
        end

        % sys: ODESys class ref
        function newLgtMdlLaTeX = updateLightAttenuationModel(sys)
            [reactorType,params] = sys.environs.(sys.activeEnv).getReactorSpecificParams();
            lgtAttnModelName = params.lgtAttnModel;
            lgtBioAttn = params.lgtBioAttn;
            lgtChemAttn = params.lgtChemAttn;

            if ~lgtBioAttn && ~lgtChemAttn
                newLgtMdlLaTeX = "I=I_{0}";
                lgtAttnModel = "I_0";
            else
                if isempty(sys.species) && isempty(sys.chemicals)
                    newLgtMdlLaTeX = "I=I_{0}";
                    lgtAttnModel = "I_0";
                else
                    compSyms = "(";
                    if lgtBioAttn
                        for k=1:1:length(sys.species)
                            if k == 1
                                compSyms = compSyms+"X_"+k;
                            else
                                compSyms = compSyms+"+X_"+k;
                            end
                        end
                    end
                    if lgtChemAttn
                        for k=1:1:length(sys.chemicals)
                            if strcmp(compSyms,"(")
                                compSyms = compSyms+"C_"+k;
                            else
                                compSyms = compSyms+"+C_"+k;
                            end                
                        end
                    end
                    compSyms = compSyms + ")";

                    switch reactorType
                        case "Vertical Continuous Stirred Tank"
                            switch lgtAttnModelName
                                case "Model 1"
                                    newLgtMdlLaTeX = "I=\frac{\int_0^L \int_0^\sqrt{L^2-y^2} I_0\,exp(-A\,(l+B\,"+compSyms+")) dxdy}{pi.*(L./2).^2}";
                                    lgtAttnModel = "integral2(@(y,x)I_0*exp(-(A*(l+B*"+compSyms+"))),0,L,0,@(y) sqrt(L.^2-y.^2))/(pi.*(L./2).^2)";
                                case "Model 2"
                            end
                        case "Horizontal Plug Flow Reactor"

                        case "Flat Panel"
                            switch lgtAttnModelName
                                case "Model 1"
                                    newLgtMdlLaTeX = "I=\frac{\int_0^L I_0\,exp(-A\,(l+B\,"+compSyms+")) dl}{L}";
                                    lgtAttnModel = "integral(@(l)I_0*exp(-(A*(l+B*"+compSyms+"))),0,L)/L";
                                case "Model 2"

                            end
                        case "Raceway Pond"
        
                        case "Shaking Flask"
        
                        case "Benchtop Bioreactor"
        
                        case "Eutrophic Lake System"


                    end
                end
            end

            sys.updateHelperFuncs("val",'Light Intensity','I',lgtAttnModel);
            sys.updateHelperFuncs("latex",'Light Intensity','I',newLgtMdlLaTeX);
            newLgtMdlLaTeX = "$" + newLgtMdlLaTeX + "$";
        end

        % sys: ODESys class ref
        function updatePBRMaxVol(sys)
            sys.PBR_max_vol = double(string(sys.environs.(sys.activeEnv).getMaxV().getSubFuncVal()));

            if sys.checkLiqMaxVol()
                % then nothing, all good
            else
                % throw an error to the user that there will be a
                % compilation or calculation error in due to the volumes
                % being outside the valid bounds
            end
        end

        % sys: ODESys class ref
        function valid = checkLiqMaxVol(sys)
            sus_solid_vol = 0;
            chems = struct2cell(sys.chemicals);
            for k=1:1:length(chems)
                if strcmp(chems{k}.getType(),'Suspended Solid Sorbent')
                    next_vol = chems{k}.getInitConc();
                    switch chems{k}.getInitConcUnit()
                        case 'L'
                            
                        case 'm^3'
                            next_vol = next_vol .* 1E3;
                        case 'mL'
                            next_vol = next_vol .* 1E-3;
                    end
                    sus_solid_vol = sus_solid_vol + next_vol;
                end
            end
            V_comp = sys.environs.(sys.activeEnv).getVComp();
            V_init = V_comp.getInitConc();
            switch V_comp.getInitConcUnit()
                case 'L'
                case 'm^3'
                    V_init = V_init .* 1E3;
                case 'mL'
                    V_init = V_init .* 1E-3;
            end

            test = sys.PBR_max_vol
            valid = sys.PBR_max_vol*0.99 >= (sus_solid_vol + V_init); % L
        end

        % sys: ODESys class ref
        function PBR_max_vol = getPBRMaxVol(sys)
            PBR_max_vol = sys.PBR_max_vol;
        end

        % sys: ODESys class ref
        function max_liq_vol = getMaxLiqVol(sys)
            max_liq_vol = sys.max_liq_vol;
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

            % removing mass transfer govFunc from P_comp
            chem = sys.chemicals.(chemName);
            P_comp = sys.environs.(sys.activeEnv).getPComp();
            P_comp.removeModel("k_"+chem.sym+"_g*"+chem.h_const_sym+"*("+chem.sym+"-"+chem.bulk_gas_sym+")", ...
                chem.sym+" liq_to_bulk",'Main');

            if strcmp(chem.getType(),'Suspended Solid Sorbent')
                helperVal = sys.getHelperFunc('Total Solvent Volume');
                helperVal = erase(helperVal,"+V_S_"+chem.getNum()); 
                sys.updateHelperFuncs("val",'Total Solvent Volume','V_tot',helperVal);

                sys.updatePBRMaxVol();
            end

            % remove HConst helpers
            if chem.is_vol
                HConstHelpers = chem.getHConstHelperFuncs();
                for k=1:1:length(HConstHelpers)
                    sys.addRmHelperFuncs(HConstHelpers{k}.funcName,HConstHelpers{k}.funcSym,HConstHelpers{k}.funcVal,false);
                end

                chem.removeModelByName('Gas MT','Liquid');
                chem.removeModelByName('Liquid MT','Gas');
            end

            % updates IOFuncs
            sys.updateIOFunc();

            sys.chemicals = rmfield(sys.chemicals, chemName);
        end

        % sys: ODESys class ref
        function [paramVals, paramUnits, funcHandles] = getCurrentEnvironParams(sys)
            paramVals = [char(string(sys.model_runtime));sys.environs.(sys.activeEnv).getParamVals(sys.getDefaultParamVals())];
            funcHandles = sys.environs.(sys.activeEnv).getParamFuncHandles(sys.getModelVarNames());
            custDefault = EnvDefaults();
            paramUnits = struct2cell(custDefault.Units);
        end

       % sys: ODESys class ref, prop: string, onlyMultiPhase: Boolean
        function res = getSoluteChems(sys,prop,onlyMultiPhase)
            if onlyMultiPhase
                res = {};
                switch prop
                    case 'name'
                        chems = struct2cell(sys.chemicals);
                        for k=1:1:length(chems)
                            % test3 = chems{k}.getName()
                            % test2 = isempty(chems{k}.sorpFuncParams)
                            if strcmp(chems{k}.type,'Chemical Solute') && ...
                                    (chems{k}.is_vol || ~isempty(chems{k}.sorpFuncParams))
                                res{k,1} = chems{k}.name; %#ok<AGROW>
                            end
                        end
                    case 'comp'
                        chems = struct2cell(sys.chemicals);
                        for k=1:1:length(chems)
                            if strcmp(chems{k}.type,'Chemical Solute') && ...
                                    (chems{k}.is_vol || isempty(chems{k}.sorpFuncParams))
                                res{k,1} = chems{k}; %#ok<AGROW>
                            end
                        end
                end
            else
                switch prop
                    case 'name'
                        chems = struct2cell(sys.chemicals);
                        res = cell(length(chems),1);
                        for k=1:1:length(chems)
                            if strcmp(chems{k}.type,'Chemical Solute') || strcmp(chems{k}.type,'Suspended Solid Sorbent')
                                res{k} = chems{k}.name;
                            end
                        end
                    case 'comp'
                        chems = struct2cell(sys.chemicals);
                        res = cell(length(chems),1);
                        for k=1:1:length(chems)
                            if strcmp(chems{k}.type,'Chemical Solute') || strcmp(chems{k}.type,'Suspended Solid Sorbent')
                                res{k} = chems{k};
                            end
                        end
                end
            end
        end

        % sys: ODESys class ref, compName: string, funcVal: string,
        % funcName: string
        function updateModel(sys, compName, funcVal, funcName)
            % get component object
            comp = sys.getCompByName(compName);

            % gets the param expressions of dependent variables
            depVars = sys.getModelVarNames(["specs","chems","envs","helpers","updateModel","matlabFuncs"]);

            % param and depVar separation for funcVal
            symVarStrArr = string(findVars(char(funcVal)));
            % creates array of system variables as strings
            varStr = string(intersect(symVarStrArr,depVars));
            % creates array of parameters as strings
            paramStr = string(setxor(symVarStrArr,varStr));

            % adding model to Comp
            % funcName = "BioGovFunc";
            comp.setModel(funcVal, funcName, paramStr, sys.getDefaultParamVals());

            for k=1:1:length(paramStr)
                sys.updateMultiParamsInv(compName,paramStr(k));
            end
        end

        % sys: ODESys class ref, compName: string, funcVal: string,
        % funcName: string, phase: string
        function addModel(sys, compName, funcVal, funcName, phase)
            % get component object
            comp = sys.getCompByName(compName);

            % gets the param expressions of dependent variables
            depVars = sys.getModelVarNames(["specs","chems","envs","helpers","updateModel","matlabFuncs"]);

            % param and depVar separation for funcVal
            symVarStrArr = string(findVars(char(funcVal)));
            % creates array of system variables as strings
            varStr = string(intersect(symVarStrArr,depVars));
            % creates array of parameters as strings
            paramStr = string(setxor(symVarStrArr,varStr));

            % adding model to Comp
            comp.addModel(funcVal, funcName, paramStr, phase, sys.getDefaultParamVals());
            
            for k=1:1:length(paramStr)
                sys.updateMultiParamsInv(compName,paramStr(k));
            end
        end

        % sys: ODESys class ref, compName: string, funcVal: string,
        % funcName: string
        function addSuspendedSolidVolumeFunc(sys,compName,funcVal,funcName)
            % get component object
            comp = sys.getCompByName(compName);

            % gets the param expressions of dependent variables
            depVars = sys.getModelVarNames(["specs","chems","envs","helpers","updateModel","matlabFuncs"]);

            % param and depVar separation for funcVal
            symVarStrArr = string(findVars(char(funcVal)));
            % creates array of system variables as strings
            varStr = string(intersect(symVarStrArr,depVars));
            % creates array of parameters as strings
            paramStr = string(setxor(symVarStrArr,varStr));

            % adding model to Comp
            comp.addSuspendedSolidVolumeFunc(funcVal, funcName, paramStr, sys.getDefaultParamVals());

            for k=1:1:length(paramStr)
                sys.updateMultiParamsInv(compName,paramStr(k));
            end
        end

        % sys: ODESys class ref, compName: string, funcVal: string, phaseA:
        % string, phaseB: string, funcType: string
        function setMTModel(sys,compName,funcVal,phaseA,phaseB,funcType)
            comp = sys.getCompByName(compName);
            
            % gets the param expressions of dependent variables
            depVars = sys.getModelVarNames(["specs","chems","envs","helpers","updateModel","matlabFuncs"]);

            % param and depVar separation for funcVal
            symVarStrArr = string(findVars(char(funcVal)));
            % creates array of system variables as strings
            varStr = string(intersect(symVarStrArr,depVars));
            % creates array of parameters as strings
            paramStr = string(setxor(symVarStrArr,varStr));

            comp.setMTModel(funcVal,phaseA,phaseB,funcType,paramStr,sys.getDefaultParamVals());

            for k=1:1:length(paramStr)
                sys.updateMultiParamsInv(compName,paramStr(k));
            end
        end

        % sys: ODESys class ref
        function defaultParamVals = getDefaultParamVals(sys,varargin)
            % should query the following parameters that aren't editable
            % within MAGMA:
            if ~isempty(varargin)
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            else
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp'), ...
                    {sys.environs.(sys.activeEnv).getH3OComp(),sys.environs.(sys.activeEnv).getOHComp()}];
            end

            % 1. MW of each component
            defaultParamVals = struct('placeholder',1);
            for k=1:1:length(comps)
                defaultParamVals.(comps{k}.MW_sym) = comps{k}.MW;
            end
            % 2. Density of each component
            try
                for k=1:1:length(comps)
                    defaultParamVals.(comps{k}.den_sym) = comps{k}.den;
                end
            catch err
                err.stack.line
            end
            % 3. Heat capacities of each component
            try
                for k=1:1:length(comps)
                    defaultParamVals.(comps{k}.Cp_sym) = comps{k}.Cp;
                    defaultParamVals.(comps{k}.Cpg_sym) = comps{k}.Cpg;
                end
            catch err
                err.stack.line
            end
        end

        % sys: ODESys class ref, prop: string
        function res = getSpecies(sys, prop)
            switch prop
                case 'name'
                    specs = struct2cell(sys.species);
                    % need to write loop for this
                    res = cell(1,length(specs));
                    for k=1:1:length(specs), res{k} = specs{k}.name; end
                case 'comp'
                    res = struct2cell(sys.species)';
            end
        end

        % sys: ODESys class ref, prop: string
        function res = getChemicals(sys, prop)
            switch prop
                case 'name'
                    chems = struct2cell(sys.chemicals);
                    % need to write loop for this
                    res = cell(1,length(chems));
                    % ext_ct = 0;
                    for k=1:1:length(chems)
                        res{k} = chems{k}.name;
                        % if chems{k}.is_vol
                        %     res{k+ext_ct+1} = [chems{k}.name,' (Gas Phase)'];
                        %     % res{i+ext_ct+2} = [chems{i+ext_ct}.name,' (Sparge Bubbles)'];
                        %     ext_ct = ext_ct + 1;
                        % end
                    end
                case 'comp'
                    res = struct2cell(sys.chemicals)';
            end
        end

        % sys: ODESys class ref
        function comps = getAllComps(sys)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
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
                sys.helperFuncs{end}.initParams(sys.getModelVarNames(),sys.getDefaultParamVals());
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
                            sys.helperFuncs{k}.initParams(sys.getModelVarNames(),sys.getDefaultParamVals());
                        end
                    end
                case "latex"
                    for k=1:1:length(sys.helperFuncs)
                        if strcmp(sys.helperFuncs{k}.getSubFuncName(),funcName) && ...
                                strcmp(sys.helperFuncs{k}.getSubFuncSym(),funcSym)
                            sys.helperFuncs{k}.setSubFuncLaTeX(funcVal);
                        end
                    end
            end
        end

        % sys: ODESys class ref, funcName: string
        function [funcSym,funcVal,funcIdx] = getHelperFunc(sys,funcName)
            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName,funcName)
                    funcSym = sys.helperFuncs{k}.getSubFuncSym();
                    funcVal = sys.helperFuncs{k}.getSubFuncVal();
                    funcIdx = k;
                    break;
                end
            end
        end

        % sys: ODESys class ref, compName: string, phase: string
        function sym = getCompSymByName(sys,compName)
            % here getting correct phase of correct comp
            % by the name that it is returned by the UI (e.g.
            % InteractCompDD)
            comp = sys.getCompByName(compName);
            sym = "";
            if any(strcmp(compName,{'Temperature','Pressure','Volume','Hydronium','Hydroxide'}))
                sym = comp.getSym();
            elseif contains(compName,' (Gas Phase)')
                sym = comp.getBulkGasSym();
            elseif contains(compName,' (Liquid Phase)')
                sym = comp.getSym();
            else
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                for k=1:1:length(comps)
                    for l=1:1:length(comps{k}.sorpFuncParams)
                        if contains(compName,comps{k}.sorpFuncParams{l}.solventName)
                            sym = comps{k}.sorpFuncParams{l}.solventSym;
                        end
                    end
                end

                if strcmp(sym,"")
                    sym = comp.getSym();
                end
            end
        end

        % sys: ODEsys class ref, compName: string, valName: string,
        % editing: boolean, editedFunc: string
        function [existFunc,existFuncUnicode] = getGovFuncByCompName(sys,compName)
            comp = sys.getCompByName(compName);

            phaseName = '';
            if any(strcmp(compName,{'Temperature','Pressure','Volume','Hydronium','Hydroxide'}))
                phase = 'Liquid Phase';
            elseif contains(compName,'(Liquid Phase)')
                phase = 'Liquid Phase';
            elseif contains(compName,'(Gas Phase)')
                phase = 'Gas Phase';
            else
                phase = 'Suspended Solid Phase';
                for k=1:1:length(comp.sorpFuncParams)
                    if contains(compName,comp.sorpFuncParams{k}.solventName)
                        phaseName = comp.sorpFuncParams{k}.solventName;
                        break;
                    end
                end
            end

            funcs = comp.getGovFuncs(phase,phaseName);
            existFunc = "";
            for k=1:1:length(funcs)
                if isempty(char(existFunc))
                    operator = "";
                else
                    operator = "+";
                end
                existFunc = existFunc + operator + funcs{k};
            end
            existFunc = uni2latex(char(existFunc));
            existFuncUnicode = funcs{1};
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

        % sys: ODESys class ref, paramSym: string,
        % newVal: number, newUnit: string, newParamName: string
        function updateMultiParams(sys,paramSym,newVal,newUnit,newParamName,convertUnits)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];

            for k=1:1:length(comps)
                fp = [comps{k}.funcParams,comps{k}.gasBulkFuncParams,comps{k}.sorpFuncParams];
                for l=1:1:length(fp)
                    for m=1:1:length(fp{l}.params)
                        if strcmp(fp{l}.params{m}.sym,paramSym)
                            funcName = fp{l}.funcName;
                            comps{k}.updateParam(paramSym,newVal,newUnit,newParamName,funcName,convertUnits);
                        end
                    end
                end
            end

            env_subfuncs = sys.environs.(sys.activeEnv).subfuncs;
            for k=1:1:length(env_subfuncs)
                envFuncName = env_subfuncs{k}.getSubFuncName();
                grthParams = sys.getGrthParamsByEnvFuncName(envFuncName);
                if isempty(grthParams), break; end
                if any(strcmp(grthParams(:,3),paramSym))
                    sys.environs.(sys.activeEnv).updateEnvFuncParams(envFuncName,paramSym,newVal,newUnit,newParamName,sys.getDefaultParamVals());
                end
            end

            for k=1:1:length(sys.helperFuncs)
                helperFuncName = sys.helperFuncs{k}.getSubFuncName();
                grthParams = sys.getGrthParamsByHelperFuncName(helperFuncName);
                if isempty(grthParams), break; end
                if any(strcmp(grthParams(:,3),paramSym))
                    sys.helperFuncs{k}.updateParams(newParamName,paramSym,newVal,newUnit,true,sys.getDefaultParamVals());
                end
            end
        end

        % sys: ODESys class ref, paramSym: string
        function updateMultiParamsInv(sys,compName,paramSym)
            comp = sys.getCompByName(compName);
            compType = comp.getType();
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];

            if any(strcmp(compType,{'Biological Solute','Chemical Solute'}))
                for k=1:1:length(comps)
                    fp = [comps{k}.funcParams,comps{k}.gasBulkFuncParams,comps{k}.sorpFuncParams];
                    for l=1:1:length(fp)
                        for m=1:1:length(fp{l}.params)
                            if strcmp(fp{l}.params{m}.sym,paramSym)
                                newVal = fp{l}.params{m}.val;
                                newParamName = fp{l}.params{m}.name;
                                newUnit = fp{l}.params{m}.unit;
                                funcName = fp{l}.funcName;
                                comp.updateParam(paramSym,newVal,newUnit,newParamName,funcName,true);
                                return;
                            end
                        end
                    end
                end
            else
                env_subfuncs = sys.environs.(sys.activeEnv).subfuncs;
                for k=1:1:length(env_subfuncs)
                    envFuncName = env_subfuncs{k}.getSubFuncName();
                    grthParams = sys.getGrthParamsByEnvFuncName(envFuncName);
                    if any(strcmp(grthParams(:,3),paramSym))
                        idx = strcmp(grthParams(:,3),paramSym);
                        newVal = grthParams(idx,5);
                        newParamName = grthParams(idx,4);
                        newUnit = grthParams(idx,6);
                        sys.environs.(sys.activeEnv).updateEnvFuncParams(envFuncName,paramSym,newVal,newUnit,newParamName,sys.getDefaultParamVals());
                    end
                end

                for k=1:1:length(sys.helperFuncs)
                    helperFuncName = sys.helperFuncs{k}.getSubFuncName();
                    grthParams = sys.getGrthParamsByHelperFuncName(helperFuncName);
                    if any(strcmp(grthParams(:,3),paramSym))
                        idx = strcmp(grthParams(:,3),paramSym);
                        newVal = grthParams(idx,5);
                        newParamName = grthParams(idx,4);
                        newUnit = grthParams(idx,6);
                        sys.helperFuncs{k}.updateParams(newParamName,paramSym,newVal,newUnit,true,sys.getDefaultParamVals());
                    end
                end
            end
        end

        % sys: ODESys class ref, prop: char, val: string
        function type = getCompType(sys, prop, val)
            % gets all species and chemical names
            specs = sys.getSpecies(prop);
            chems = sys.getChemicals(prop);

            if val == "Temperature"
                type = "temp";
            elseif val == "Pressure"
                type = "press";
            elseif val == "Volume"
                type = "vol";
            elseif val == "Hydronium"
                type = "acid";
            elseif val == "Hydroxide"
                type = "base";
            elseif any(contains(val,specs))
                type = "Biological Solute";
            elseif any(contains(val,chems))
                type = "Chemical";
            else
                type = "helper";
            end
        end

        % sys: ODESys class ref, compName: string
        function comp = getCompByName(sys,compName)
            compType = sys.getCompType('name',compName);
            comp = [];
            if compType == "Biological Solute"
                specs = sys.getSpecies('comp');
                for k=1:1:length(specs)
                    if contains(compName,specs{k}.getName())
                        comp = specs{k};
                    end
                end
            elseif contains(compType,"Chemical")
                chems = sys.getChemicals('comp');
                for k=1:1:length(chems)
                    idx = regexp(compName,'\(');
                    compName = char(compName);
                    if ~isempty(idx), compName = compName(1:idx-2); end
                    if contains(compName,chems{k}.getName())
                        comp = chems{k};
                    end
                end
            elseif compType == "temp"
                comp = sys.environs.(sys.activeEnv).getTComp();
            elseif compType == "press"
                comp = sys.environs.(sys.activeEnv).getPComp();
            elseif compType == "vol"
                comp = sys.environs.(sys.activeEnv).getVComp();
            elseif compType == "acid"
                comp = sys.environs.(sys.activeEnv).getH3OComp();
            elseif compType == "base"
                comp = sys.environs.(sys.activeEnv).getOHComp();
            elseif compType == "helper"
                for k=1:1:length(sys.helperFuncs)
                    if strcmp(compName,sys.helperFuncs{k}.getSubFuncName())
                        comp = sys.helperFuncs{k};
                    end
                end
            else
                comp = [];
                disp("Doesn't Work, ODESys.m: line 605");
                return;
            end
        end

        % sys: ODESys class ref, type: string
        function names = getCompNamesByType(sys,type,allPhases)
            if allPhases
                if strcmp(type,"Biological Solute")
                    specs = sys.getSpecies('comp');
                    names = {};
                    for k=1:1:length(specs)
                        names{end+1} = [char(specs{k}.getName()),' (Liquid Phase)']; %#ok<AGROW>
                    end
                elseif strcmp(type,"Chemical Solute")
                    chems = sys.getChemicals('comp');
                    names = {};
                    for k=1:1:length(chems)
                        names{end+1} = [char(chems{k}.getName()),' (Liquid Phase)']; %#ok<AGROW>
                        if chems{k}.is_vol
                            names{end+1} = [char(chems{k}.getName()),' (Gas Phase)']; %#ok<AGROW>
                        end
                        for l=1:1:length(chems{k}.sorpFuncParams)
                            names{end+1} = [char(chems{k}.getName()),' (',chems{k}.sorpFuncParams{l}.solventName,')']; %#ok<AGROW>
                        end
                    end
                elseif strcmp(type,"Environ")
                    names = {};
                    comps = sys.environs.(sys.activeEnv).getAllEnvComps();
                    for k=1:1:length(comps)
                        names{k} = comps{k}.getName(); %#ok<AGROW>
                    end
                end
            else
                if strcmp(type,"Biological Solute")
                    specs = sys.getSpecies('comp');
                    names = {};
                    for k=1:1:length(specs)
                        names{end+1} = char(specs{k}.getName()); %#ok<AGROW>
                    end
                elseif strcmp(type,"Chemical Solute")
                    chems = sys.getChemicals('comp');
                    names = {};
                    for k=1:1:length(chems)
                        names{end+1} = char(chems{k}.getName()); %#ok<AGROW>
                    end
                elseif strcmp(type,"Environ")
                    names = {};
                    comps = sys.environs.(sys.activeEnv).getAllEnvComps();
                    for k=1:1:length(comps)
                        names{k} = comps{k}.getName(); %#ok<AGROW>
                    end
                end
            end
        end

        % sys: ODESys class ref
        function names = getEnvFuncNames(sys)
            names = sys.environs.(sys.activeEnv).getParamNames();
        end

        % sys: ODESys class ref
        function [reactorType,params] = getReactorSpecificParams(sys)
            [reactorType,params] = sys.environs.(sys.activeEnv).getReactorSpecificParams();
        end

        % sys: ODESys class ref
        function names = getHelperFuncNames(sys)
            names = string.empty(0,1);
            for k=1:1:length(sys.helperFuncs)
                names(k,1) = sys.helperFuncs{k}.getSubFuncName();
            end
        end

        % sys: ODESys class ref, name: string
        function helperSym = getHelperSymByName(sys,name)
            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName(),name)
                    helperSym = sys.helperFuncs{k}.getSubFuncSym();
                end
            end
        end

        % sys: ODESys class ref, name: string
        function helperFuncVal = getHelperFuncValByName(sys,name)
            for k=1:1:length(sys.helperFuncs)
                if strcmp(sys.helperFuncs{k}.getSubFuncName(),name)
                    helperFuncVal = uni2latex(char(sys.helperFuncs{k}.getSubFuncVal()));
                end
            end
        end

        % sys: ODESys class ref
        function names = getEnvAndHelperFuncNames(sys)
            names = string.empty(0,1);
            helpers = [sys.helperFuncs,sys.environs.(sys.activeEnv).subfuncs];
            for k=1:1:length(helpers)
                names(k,1) = helpers{k}.getSubFuncName();
            end
        end

        % sys: ODESys class ref, name: string
        function helperSym = getEnvAndHelperSymByName(sys,name)
            helpers = [sys.helperFuncs,sys.environs.(sys.activeEnv).subfuncs];
            for k=1:1:length(helpers)
                if strcmp(helpers{k}.getSubFuncName(),name)
                    helperSym = helpers{k}.getSubFuncSym();
                end
            end
        end

        % sys: ODESys class ref, name: string
        function helperFuncVal = getEnvAndHelperFuncValByName(sys,name)
            helpers = [sys.helperFuncs,sys.environs.(sys.activeEnv).subfuncs];
            for k=1:1:length(helpers)
                if strcmp(helpers{k}.getSubFuncName(),name)
                    if strcmp(helpers{k}.getSubFuncName(),'Light Intensity')
                        helperFuncVal = helpers{k}.getSubFuncLaTeX();
                    else
                        helperFuncVal = uni2latex(char(helpers{k}.getSubFuncVal()));
                    end
                end
            end
        end

        % sys: ODESys class ref
        % Model Param Nomenclature:
        % Species concentrations: X1, X2, ... , Xn
        % Chemical concentrations: C1, C2, ... , Cn
        % Light: I
        % Time: t
        % Temperature: T
        % Pressure: P
        % Environment Volume: V
        % this method will return the model var names of all biological and
        % chemical solute species in all phases
        % intended for CompParamDD.Items and SysVarT.Data,
        % RegCompParamDD.Items, and RegSysVarDD.Items
        function vars = getModelVarNames(sys,varargin)
            % compiling all parameters from Environment and Species and
            % Chemical concentrations to push into SysVarT
            specs = struct2cell(sys.species);
            chems = struct2cell(sys.chemicals);

            if isempty(varargin)
                names_req = ["specs","chems","envs","helpers","updateModel"];
            elseif strcmp(varargin{1},"plot")
                names_req = "plot";
            else
                names_req = varargin{1};
            end

            vars = {};
            if any(strcmp(names_req,"envs"))
                % Environmental parameters
                envParams = Environment.getParamNames();
                for k=1:1:size(envParams,1)
                    vars{end+1,1} = envParams{k,1}; %#ok<AGROW>
                    vars{end,2} = envParams{k,2};
                end
            end

            if any(strcmp(names_req,"specs"))
                vars = sys.getModelSpecVarNames(specs,vars);
            end

            if any(strcmp(names_req,"chems"))
                vars = sys.getModelChemVarNames(chems,vars);
            end

            if any(strcmp(names_req,"helpers"))
                % Subfuncs
                for k=1:1:length(sys.helperFuncs)
                    vars{end+1,1} = char(sys.helperFuncs{k}.getSubFuncName()); %#ok<AGROW>
                    vars{end,2} = sys.helperFuncs{k}.getSubFuncSym();
                end
            end

            if any(strcmp(names_req,"defaultParams"))
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp'), ...
                    {sys.environs.(sys.activeEnv).getH3OComp(),sys.environs.(sys.activeEnv).getOHComp()}];
                % uneditable parameters
                for k=1:1:length(comps)
                    vars{end+1,1} = [comps{k}.getName(),' Molecular Weight']; %#ok<AGROW>
                    vars{end,2} = comps{k}.MW_sym;
                end
            end

            if any(strcmp(names_req,"updateModel"))
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
                for k=1:1:length(comps)
                    if strcmp(comps{k}.getType(),'Suspended Solid Sorbent')
                        vars{end+1,1} = [char(comps{k}.getSym()),' Inlet Volumetric Flow Rate']; %#ok<AGROW>
                        vars{end,2} = ['L_s_i_',char(comps{k}.sym)];
                        vars{end+1,1} = [char(comps{k}.getSym()),' Outlet Volumetric Flow Rate']; %#ok<AGROW>
                        vars{end,2} = ['L_s_o_',char(comps{k}.sym)];
                    else
                        vars{end+1,1} = [char(comps{k}.getSym()),' Mass Flow In']; %#ok<AGROW>
                        vars{end,2} = ['D_i_',char(comps{k}.sym)];
                        vars{end+1,1} = [char(comps{k}.getSym()),' Liquid Concentration In']; %#ok<AGROW>
                        vars{end,2} = [char(comps{k}.getSym()),'_i'];
                        if comps{k}.is_vol
                            vars{end+1,1} = [char(comps{k}.getSym()),' Inlet Bubble Gas Mol Fraction']; %#ok<AGROW>
                            vars{end,2} = ['Y_',char(comps{k}.getSym()),'_i'];
                        end
                        for l=1:1:length(comps{k}.sorpFuncParams)
                            vars{end+1,1} = [char(comps{k}.sorpFuncParams{l}.solventSym),' Phase Inlet Concentration']; %#ok<AGROW>
                            vars{end,2} = [char(comps{k}.sorpFuncParams{l}.solventSym),'_i'];
                        end
                    end
                end
                % I/O
                vars{end+1,1} = 'Outlet Pressure Change';
                vars{end,2} = 'dP_o';
                vars{end+1,1} = 'Liquid Mass Flow In';
                vars{end,2} = 'L_i';
                vars{end+1,1} = 'Liquid Mass Flow Out';
                vars{end,2} = 'L_o';
                vars{end+1,1} = 'Gas Mass Flow In';
                vars{end,2} = 'D_i';
                % ### STARTHERE: need to fix this so each suspended solid
                % phase has its own liquid in and out flow
                % vars{end+1,1} = 'Liquid Mass Flow In';
                % vars{end,2} = 'L_s_i';
                % vars{end+1,1} = 'Liquid Mass Flow Out';
                % vars{end,2} = 'L_s_o';
            end

            if any(strcmp(names_req,"matlabFuncs"))
                vars{end+1,1} = 'Log Base 10';
                vars{end,2} = 'log10';
                vars{end+1,1} = 'Natural Log';
                vars{end,2} = 'ln';
                vars{end+1,1} = 'Natural Exponent';
                vars{end,2} = 'exp';
            end

            if any(strcmp(names_req,"plot"))
                vars{end+1,1} = 'Model Runtime';
                vars{end,2} = 't';

                vars = sys.getModelSpecVarNames(specs,vars);
                vars = sys.getModelChemVarNames(chems,vars);

                % Environmental parameters
                envParams = Environment.getParamNames();
                for k=1:1:size(envParams,1)
                    vars{end+1,1} = envParams{k,1}; %#ok<AGROW>
                    vars{end,2} = envParams{k,2};
                end

                % Subfuncs
                for k=1:1:length(sys.helperFuncs)
                    vars{end+1,1} = char(sys.helperFuncs{k}.getSubFuncName()); %#ok<AGROW>
                    vars{end,2} = sys.helperFuncs{k}.getSubFuncSym();
                end
            end

            if any(strcmp(names_req,"nonNeg"))
                vars = sys.getModelSpecVarNames(specs,vars);
                vars = sys.getModelChemVarNames(chems,vars);

                for k=1:1:length(chems)
                    if strcmp(chems{k}.getType(),'Suspended Solid Sorbent')
                        vars{end+1,1} = [char(chems{k}.getName()),' Volume']; %#ok<AGROW>
                        vars{end,2} = chems{k}.getSym();
                    end
                end

                % Environmental parameters
                envParams = Environment.getParamNames();
                for k=1:1:size(envParams,1)
                    vars{end+1,1} = envParams{k,1}; %#ok<AGROW>
                    vars{end,2} = envParams{k,2};
                end
            end
        end

        % sys: ODESys class ref
        function vars = getModelVarNamesReg(sys)
            specs = struct2cell(sys.species);
            chems = struct2cell(sys.chemicals);

            vars = {};
            vars = sys.getModelSpecVarNames(specs,vars);

            vars = sys.getModelChemVarNames(chems,vars);

            % Environmental parameters
            % envParams = Environment.getParamNames();
            env_comps = sys.environs.(sys.activeEnv).getAllEnvComps();
            for k=1:1:length(env_comps)
                vars{end+1,1} = env_comps{k}.getName(); %#ok<AGROW>
                vars{end,2} = env_comps{k}.getSym();
            end
        end

        % sys: ODESys class ref, specs: {}, vars: {}
        function vars = getModelSpecVarNames(~,specs,vars)
            % Species Concentrations
            if length(specs) > 0 %#ok<ISMT> 
                for k=1:1:length(specs)
                    vars{end+1,1} = [char(specs{k}.name),' (Liquid Phase)']; %#ok<AGROW>
                    vars{end,2} = char(specs{k}.getSym());
                end
            end
        end

        % sys: ODESys class ref, chems: {}, vars: {}
        function vars = getModelChemVarNames(sys,chems,vars)
            % Chemical Concentrations
            if length(chems) > 0 %#ok<ISMT> 
                for k=1:1:length(chems)
                    if strcmp(chems{k}.getType(),'Chemical Solute')
                        vars{end+1,1} = [chems{k}.name,' (Liquid Phase)']; %#ok<AGROW>
                        vars{end,2} = char(chems{k}.getSym());
                        if chems{k}.is_vol
                            vars{end+1,1} = [chems{k}.name,' (Gas Phase)']; %#ok<AGROW>
                            vars{end,2} = char(chems{k}.bulk_gas_sym);
                            if sys.has_gas_in && ~any(strcmp(names_req,"plot"))
                                vars{end+1,1} = [chems{k}.name,' (Inlet Bubble Gas Mol Fraction)']; %#ok<AGROW>
                                vars{end,2} = ['Y_i_',char(chems{k}.sym)];
                            end
                        end
                        solvents = {};
                        for l=1:1:length(chems{k}.sorpFuncParams)
                            if ~any(strcmp(solvents,chems{k}.sorpFuncParams{l}.solventName))
                                vars{end+1,1} = [chems{k}.name,' (',chems{k}.sorpFuncParams{l}.solventName,' Phase)']; %#ok<AGROW>
                                vars{end,2} = chems{k}.sorpFuncParams{l}.solventSym;
                                solvents{end+1} = chems{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end

        % sys: ODESys class ref, name: string
        function params = getGrthParamsByFuncName(sys,funcName)
            subfuncs = [sys.environs.(sys.activeEnv).subfuncs,sys.helperFuncs];
            subfuncNames = string(zeros(size(subfuncs)));
            for k=1:1:length(subfuncs), subfuncNames(k) = subfuncs{k}.getSubFuncName(); end
            if any(strcmp(subfuncNames,funcName))
                params = sys.getGrthParamsByHelperFuncName(funcName);
            else
                phase = 'Main';
                comp = sys.getCompByName(funcName);
                if contains(funcName,'(Gas Phase)')
                    phase = 'Gas';
                elseif contains(funcName,'(Liquid Phase)')
                    phase = 'Liquid';
                else
                    comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                    for k=1:1:length(comps)
                        for l=1:1:length(comps{k}.sorpFuncParams)
                            if contains(funcName,comps{k}.sorpFuncParams{l}.solventName)
                                phase = comps{k}.sorpFuncParams{l}.solventName;
                            end
                        end
                    end
                end
                params = comp.getGrthParams(phase);
            end
        end

        % sys: ODESys class ref, envFuncName: string
        function params = getGrthParamsByEnvFuncName(sys,envFuncName)
            switch envFuncName
                case "Temperature"
                    params = sys.environs.(sys.activeEnv).getTComp().getGrthParams('Main');
                case "Pressure"
                    params = sys.environs.(sys.activeEnv).getPComp().getGrthParams('Main');
                case "Volume"
                    params = sys.environs.(sys.activeEnv).getVComp().getGrthParams('Main');
                case "Incident Light"
                    params = sys.environs.(sys.activeEnv).getGrthParamsByEnvFuncName(envFuncName);
                case "Maximum Volume"
                    params = sys.environs.(sys.activeEnv).getGrthParamsByEnvFuncName(envFuncName);
                case "Hydronium"
                    params = sys.environs.(sys.activeEnv).getH3OComp().getGrthParams('Main');
                case "Hydroxide"
                    params = sys.environs.(sys.activeEnv).getOHComp().getGrthParams('Main');
            end
        end

        % sys: ODESys class ref, helperFuncName: string
        function params = getGrthParamsByHelperFuncName(sys,helperFuncName)
            subfuncs = [sys.environs.(sys.activeEnv).subfuncs,sys.helperFuncs];
            for k=1:1:length(subfuncs)
                if strcmp(subfuncs{k}.getSubFuncName(),helperFuncName)
                    paramNames = subfuncs{k}.getSubFuncParamNames();
                    paramNums = 1:1:length(paramNames);
                    if isempty(paramNums), params = {}; return; end
                    paramSyms = subfuncs{k}.getSubFuncParamSyms();
                    paramVals = subfuncs{k}.getSubFuncParamVals();
                    paramUnits = subfuncs{k}.getSubFuncParamUnits();
                    params = cell(length(paramNums),6);
                    for l=paramNums
                        val = unit_standardization(paramVals(l),paramUnits(l));
                        params{l,1} = char(string(paramNums(l)));
                        params{l,2} = helperFuncName;
                        params{l,3} = char(paramSyms(l));
                        params{l,4} = char(paramNames(l));
                        params{l,5} = val;
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

        % sys: ODESys class ref, compName: string, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, lowLimOp: string,
        % upLimOp: string, phaseA: string, phaseB: string
        function pwTInfo = addNewPWExpr(sys,compName,newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB)
            comp = sys.getCompByName(compName);
            pwTInfo = comp.addNewPWExpr(newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB);
        end

        % sys: ODESys class ref, compName: string, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, lowLimOp: string,
        % upLimOp: string, phaseA: string, phaseB: string
        function pwTInfo = editPWExpr(sys,compName,newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB)
            comp = sys.getCompByName(compName);
            pwTInfo = comp.editPWExpr(newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB); 
        end

        % sys: ODESys class ref, compName: string, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, lowLimOp: string,
        % upLimOp: string, phaseA: string, phaseB: string
        function pwTInfo = rmPWExpr(sys,compName,newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB)
            comp = sys.getCompByName(compName);
            pwTInfo = comp.rmPWExpr(newExpr,lowLim,upLim,sysVar,elseVar,lowLimOp,upLimOp,phaseA,phaseB);
        end


        % ### FIXME: plot specifications
        % need to add code to allow user to specify a t range and precision
        % of t

        % sys: ODEsys class ref, tRange: number[], tPtsNb: number
        % ### ACTUALLYSTARTHERE: if volume of gas phase is 0, it causes
        % errors in liquid gas mass transfer equations; how to resolve
        % this?
        % Handle this by:
        %   if V_liq+sum(V_sus_solids) >= 0.99*V_max: shut off all input streams,
        %   show warning at end of run
        %   if V_liq <= 0.01*V_max: shut off all liquid output streams, show warning at end
        %   of run
        % Also remove Volume and Pressure rate equation modification from
        % GovFuncs page
        function compileModel(sys)
            % 1. loop through funcParams in each Component and assemble
            % functions into full governing function
            %   replace each parameter name with p({gloNum})
            %   can get the parameters to be converted from string to sym
            %   with all parameters saved in the p array by:
            %   also need to convert names of system variables
            % 2. convert each governing function to sym function
            % 3. create ODESys of sym functions

            % ### FIXME: need to test to make sure it works
            % ### UPDATE: break into pieces so MAGMA can do one piece at a
            % time as needed
            % ### STARTHERE: add functionality to compile I/O stream flows
            % into model: need to create code to generate v() array

            try
                sys.dydt = "@(t,y,p,f,v) [";
                sys.param = [];
                sys.f = {};
                sys.v = [];
    
                comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                sys.var_out_key = cell(length(comps),2);
                sys.S_nums = {};
                ext_ct = 0;
                for k=1:1:length(comps)
                    sys.var_out_key{k+ext_ct,1} = comps{k}.getName();
                    sys.var_out_key{k+ext_ct,2} = comps{k}.getSym();
                    if comps{k}.is_vol  
                        ext_ct = ext_ct + 1;
                        sys.var_out_key{k+ext_ct,2} = comps{k}.getName();
                        sys.var_out_key{k+ext_ct,2} = comps{k}.getBulkGasSym();
                    end
                    for l=1:1:length(comps{k}.sorpFuncParams)
                        ext_ct = ext_ct + 1;
                        sys.var_out_key{k+ext_ct,2} = comps{k}.getName();
                        sys.var_out_key{k+ext_ct,2} = comps{k}.sorpFuncParams{l}.solventSym;
                    end
                    if strcmp(comps{k}.getType(),'Suspended Solid Sorbent')
                        sys.S_nums{end+1} = comps{k}.getNum();
                    end
                end
    
                in_syms = sys.getModelVarNames("updateModel");
                in_syms = in_syms(:,2);
                sys.var_in_key = {};
                num_spec = length(sys.getSpecies('comp'));
                num_spec_chem = length([sys.getSpecies('comp'),sys.getChemicals('comp')]);
    
                % ### UPDATE: updating V_m value, will later find better way to
                % do this: either use a helper function, or automatically
                % incorporate this value into a function
                % sys.updateMultiParams('V_m',sys.PBR_max_vol,'L','Maximum Volume');
                subfuncs = [sys.environs.(sys.activeEnv).subfuncs,sys.helperFuncs];
                for k=1:1:length(comps) % should be length(comps) + length(environSubfuncs)
                    % add logic to compile component gov funcs in order that
                    % allows for growth-associated product funcs to be a
                    % function of biomass funcs, and substrate funcs to be a
                    % function of biomass and product funcs
                    % ### FIXME: requires special logic/restrictions in the
                    % governing function portion
    
                    % if component is_vol == true:
                    % in comp.compileGovFunc() do a pre-compilation of
                    % all of the components, then keep the compileModel portion
                    % the same
                    [govFunc, params, ~] = comps{k}.compileGovFunc(length(sys.param),sys.reg_param_ct,{},false);
                    for l=1:1:length(subfuncs)
                        if ~strcmp(subfuncs{l}.getSubFuncSym(),"t")
                            govFunc = regexprep(govFunc,subfuncs{l}.getSubFuncSym(),"f{"+l+"}(t,y,p,f,v)");
                        end
                    end
                    % loop through all in_syms and replace with v(#), and track
                    % which # is for which syms
                    for l=1:1:length(in_syms)
                        if ~isempty(regexp(govFunc,in_syms{l},'once'))
                            sys.var_in_key{end+1} = in_syms{l};
                            govFunc = replace(govFunc,in_syms{l},['v(',char(string(length(sys.var_in_key))),')']);
                        end
                    end
                    ext_ct = 0;
                    for l=1:1:length(comps)
                        if l <= num_spec
                            govFunc = regexprep(govFunc,"X_"+string(comps{l}.getNum()),"y("+l+")");
                        elseif l > num_spec && l <= num_spec_chem
                            govFunc = replace(govFunc,['C_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_C_',char(string(comps{l}.getNum()))]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['C_',char(string(comps{l}.getNum())),'_']);
                            if comps{l}.is_vol
                                ext_ct = ext_ct + 1;
                                govFunc = regexprep(govFunc,"P_"+(comps{l}.getNum()),"y("+string(l+ext_ct)+")");
                            end
                            phases = comps{l}.getSorpPhases();
                            for m=1:1:length(phases)
                                ext_ct = ext_ct + 1;
                                govFunc = regexprep(govFunc,phases{m},"y("+(l+ext_ct)+")");
                                govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_',phases{m}]);
                                govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],[phases{m},'_']);
                            end
                            if strcmp(comps{l}.getType(),'Suspended Solid Sorbent')
                                govFunc = replace(govFunc,['V_S_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                                govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_V_S_',char(string(comps{l}.getNum()))]);
                                govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['V_S_',char(string(comps{l}.getNum())),'_']);
                            end
                        elseif l == num_spec_chem + 1
                            govFunc = regexprep(govFunc,'T',['y(',char(string(l+ext_ct)),')']);
                        elseif l == num_spec_chem + 2
                            govFunc = regexprep(govFunc,'P',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'P_');
                        elseif l == num_spec_chem + 3
                            govFunc = regexprep(govFunc,'V',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'V_');
                        elseif l == num_spec_chem + 4
                            govFunc = replace(govFunc,'H3O',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_H3O');
                        elseif l == num_spec_chem + 5
                            govFunc = replace(govFunc,'OH',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_OH');
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
                for k=1:1:length(subfuncs)
                    ext_ct = 0;
                    subfuncParamCt = 1;
                    [govFunc,params,~] = subfuncs{k}.compileSubFunc(length(sys.param),sys.reg_param_ct,{},false);
                    for l=1:1:length(comps)
                        if l <= num_spec
                            govFunc = replace(govFunc,"X_"+string(comps{l}.getNum()),"y("+l+")");
                        elseif l > num_spec && l <= num_spec_chem
                            govFunc = replace(govFunc,['C_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_C_',char(string(comps{l}.getNum()))]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['C_',char(string(comps{l}.getNum())),'_']);
                            if comps{l}.is_vol
                                ext_ct = ext_ct + 1;
                                govFunc = replace(govFunc,"P_"+(comps{l}.getNum()),"y("+string(l+ext_ct)+")");
                            end
                            phases = comps{l}.getSorpPhases();
                            for m=1:1:length(phases)
                                ext_ct = ext_ct + 1;
                                govFunc = replace(govFunc,phases{m},"y("+(l+ext_ct)+")");
                                govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_',phases{m}]);
                                govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],[phases{m},'_']);
                            end
                            if strcmp(comps{l}.getType(),'Suspended Solid Sorbent')
                                govFunc = replace(govFunc,['V_S_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                                govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_V_S_',char(string(comps{l}.getNum()))]);
                                govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['V_S_',char(string(comps{l}.getNum())),'_']);
                                sys.IO_flowrates.(['V_S_',char(string(comps{l}.getNum())),'_idx']) = (l+ext_ct);
                            end
                        elseif l == num_spec_chem + 1
                            govFunc = regexprep(govFunc,'T',['y(',char(string(l+ext_ct)),')']);
                        elseif l == num_spec_chem + 2
                            govFunc = regexprep(govFunc,'P',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'P_');
                        elseif l == num_spec_chem + 3
                            govFunc = regexprep(govFunc,'V',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'V_');
                            sys.IO_flowrates.V_idx = (l+ext_ct);
                        elseif l == num_spec_chem + 4
                            govFunc = replace(govFunc,'H3O',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_H3O');
                        elseif l == num_spec_chem + 5
                            govFunc = replace(govFunc,'OH',['y(',char(string(l+ext_ct)),')']);
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_OH');
                        end
                    end
                    for l=1:1:length(subfuncs)
                        if ~strcmp(subfuncs{l}.getSubFuncSym(),"t")
                            govFunc = regexprep(govFunc,subfuncs{l}.getSubFuncSym(),"f{"+l+"}(t,y,p,f,v)");
                        end
                    end
                    % making sure the in_syms are tracked accurately from
                    % component govFuncs
                    for l=1:1:length(in_syms)
                        if ~isempty(regexp(govFunc,in_syms{l},'once'))
                            if any(strcmp(sys.var_in_key,in_syms{l}))
                                govFunc = replace(govFunc,in_syms{l},['v(',char(string(find([sys.var_in_key{:}] == in_syms{l}))),')']);
                            else
                                sys.var_in_key{end+1} = in_syms{l};
                                govFunc = replace(govFunc,in_syms{l},['v(',char(string(length(sys.var_in_key))),')']);
                            end
                        end
                    end
    
                    syms = subfuncs{k}.getSubFuncParamSyms();
                    govFuncLength = strlength(govFunc);
                    for l=1:1:length(syms)
                        for m=1:1:length(govFunc)
                            if strlength(govFunc(m)) > 1 || govFuncLength == 1
                                govFunc(m) = regexprep(govFunc(m),syms(l),"#("+(length(sys.param)+subfuncParamCt)+")");
                            end
                        end
                        subfuncParamCt = subfuncParamCt + 1;
                    end
                    govFunc = replace(govFunc,'#','p');
                    funcArgText = "@(t,y,p,f,v)";
                    sys.f{end+1,1} = subfuncs{k}.getSubFuncName();
                    sys.f{end,2} = subfuncs{k}.getSubFuncSym();
                    % test = funcArgText+govFunc
                    sys.f{end,3} = str2func(funcArgText+govFunc);
                    sys.param = [sys.param,params];
                end

                % ### STARTHERE: calculate values for v() array
                comps_v = [sys.getSpecies('comp'),sys.getChemicals('comp')];
                for k=1:1:length(sys.var_in_key)
                    switch sys.var_in_key{k}
                        case 'dP_o'
                            

                        case 'L_i'
                            sys.IO_flowrates.L_i = 0;

                            sys.IO_flowrates.L_i = sys.v(k);
                            sys.IO_flowrates.L_i_idx = k;
                        case 'L_o'
                            sys.IO_flowrates.L_o = 0;
                            for l=1:1:length(sys.output_streams)
                                [~,~,~,solventName,flowrate,flowrateU] = sys.output_streams{l}.getCompsData();
                                if strcmp(sys.output_streams{l}.getPhase(),'L') && strcmp(solventName,'Water (H2O)')
                                    switch flowrateU
                                        case 'm^3/s'
                                            flowrate = flowrate_sum .* 1000;
                                        case 'L/s'
                                        case 'mL/s'
                                            flowrate = flowrate_sum ./ 1000;
    
                                    end
                                    sys.v(k) = sys.v(k) + flowrate;
                                end
                            end
                            sys.IO_flowrates.L_o = sys.v(k);
                            sys.IO_flowrates.L_o_idx = k;
                        case 'D_i'
                            for l=1:1:length(sys.input_streams)
                                if strcmp(sys.input_streams{l}.getPhase(),'G')
                                    [compsData,basis,basisU,~,flowrate,flowrateU] = sys.input_streams{l}.getCompsData();
                                    if strcmp(basis,'flow')
                                        flowrate_sum = sum(double(string(compsData(:,2))));
                                        switch basisU
                                            case 'MPa/s'
                                                flowrate = flowrate_sum .* 1E6;
                                            case 'kPa/s'
                                                flowrate = flowrate_sum .* 1000;
                                            case 'Pa/s'
    
                                        end
                                    elseif strcmp(basis,'frac')
                                        switch flowrateU
                                            case 'MPa/s'
                                                flowrate = flowrate_sum .* 1E6;
                                            case 'kPa/s'
                                                flowrate = flowrate_sum .* 1000;
                                            case 'Pa/s'
    
                                        end
                                    end
                                    sys.v(k) = sys.v(k) + flowrate;
                                end
                            end
                        otherwise
                            if contains(sys.var_in_key{k},'L_s_i_')
                                sys.IO_flowrates.(sys.var_in_key{k}) = 0;

                                sys.IO_flowrates.(sys.var_in_key{k}) = sys.v(k);
                                sys.IO_flowrates.(string(sys.var_in_key{k})+"_idx") = k;                                
                            elseif contains(sys.var_in_key{k},'L_s_o_')
                                sys.IO_flowrates.(sys.var_in_key{k}) = 0;

                                sys.IO_flowrates.(sys.var_in_key{k}) = sys.v(k);
                                sys.IO_flowrates.(string(sys.var_in_key{k})+"_idx") = k;
                            elseif contains(sys.var_in_key{k},'D_i_')

                            elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'C_') && ~contains(sys.var_in_key{k},'D_') && ~contains(sys.var_in_key{k},'L_s_')

                            elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'Y_')

                            elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'S_') && ~contains(sys.var_in_key{k},'D_') && ~contains(sys.var_in_key{k},'L_s_')

                            end
                    end
                end
    
                % converting to function_handle
                sys.dydt = str2func(sys.dydt');
                test7 = sys.dydt
                
                % setting initial conditions
                sys.setInitCond();

                % setting solver events
                % sys.setSolverEvents();
            catch err
                err.stack.line
            end
        end

        % set initial conditions for all phases of
        % all components, evaluate the equilibrium helper functions of each
        % component
        % sys: ODESys class ref
        function setInitCond(sys)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
            sys.y0 = zeros(1,length(comps));
            sys.solverNonNeg = zeros(1,length(comps));
            ext_ct = 0;
            solvent_total_V = 0;
            for k=1:1:(length(comps))
                sys.y0(k+ext_ct) = double(string(comps{k}.getInitConc()));
                if any(strcmp(sys.nonNegSysVars,[char(comps{k}.getName()),' (Liquid Phase)'])), sys.solverNonNeg(k+ext_ct) = 1; end
                if strcmp(comps{k}.getType(),'vol') || strcmp(comps{k}.getType,'Suspended Solid Sorbent')
                    solvent_total_V = solvent_total_V + comps{k}.getInitConc();
                end
                if comps{k}.is_vol
                    ext_ct = ext_ct + 1;
                    helperName = ['Gas Phase Partial Pressure of ',char(comps{k}.getName()),' in Equilibrium with Liquid Phase'];
                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                    % ### FIXME?: unable to access initial conditions for
                    % variables inputted from Simulink, will assume 0 for
                    % all
                    sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                    if any(strcmp(sys.nonNegSysVars,[char(comps{k}.getName()),' (Gas Phase)'])), sys.solverNonNeg(k+ext_ct) = 1; end
                end
                solvents = {};
                if ~isempty(comps{k}.sorpFuncParams)
                    for l=1:1:length(comps{k}.sorpFuncParams)
                        if ~any(strcmp(solvents,comps{k}.sorpFuncParams{l}.solventName))
                            ext_ct = ext_ct + 1;
                            helperName = [char(comps{k}.name),' Concentration in ',comps{k}.sorpFuncParams{l}.solventName,' Phase in Equilibrium with Liquid Phase'];
                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                            sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                            if any(strcmp(sys.nonNegSysVars,[comps{k}.getName(),' (',char(comps{k}.sorpFuncParams{l}.solventName),' Phase)'])), sys.solverNonNeg(k+ext_ct) = 1; end
                            solvents{end+1} = comps{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                        end
                    end
                end
            end
            if (solvent_total_V + 0.01*sys.PBR_max_vol) > sys.PBR_max_vol
                % display error on volumes to users
                return;
            end
        end

        % sys: ODESys class ref, sysVar: string, add: boolean
        function updateNonNegSysVarList(sys,sysVar,add)
            if add
                sys.nonNegSysVars{end+1} = sysVar;
            else
                for k=1:1:length(sys.nonNegSysVars)
                    if strcmp(sys.nonNegSysVars{k},sysVar)
                        sys.nonNegSysVars(k) = [];
                        return;
                    end
                end
            end
        end

        % sys: ODESys class ref, solverName: string
        function updateODESolverName(sys,solverName)
            sys.solver.Solver = solverName;
        end

        % sys: ODESys class ref, solverName: string, solverSpecs: struct
        function updateODESolverSpecs(sys,solverName,solverSpecs)
            if solverName == "auto" || solverName == "stiff" || solverName == "nonstiff"
                return;
            elseif solverName == "ode23" || solverName == "ode45" || solverName == "ode78" || solverName == "ode89" || solverName == "ode113"
                sys.solver.InitialStep = solverSpecs.InitialStep;
                sys.solver.MaxStep = solverSpecs.MaxStep;
                sys.solver.NormControl = solverSpecs.NormControl;
            elseif solverName == "ode15s"
                sys.solver.InitialStep = solverSpecs.InitialStep;
                sys.solver.MaxStep = solverSpecs.MaxStep;
                sys.solver.NormControl = solverSpecs.NormControl;
                sys.solver.Vectorized = solverSpecs.Vectorized;
                sys.solver.BDF = solverSpecs.BDF;
                sys.solver.MaxOrder = solverSpecs.MaxOrder;
            elseif solverName == "ode23s" || solverName == "odet" || solverName == "ode23tb"
                sys.solver.InitialStep = solverSpecs.InitialStep;
                sys.solver.MaxStep = solverSpecs.MaxStep;
                sys.solver.NormControl = solverSpecs.NormControl;
                sys.solver.Vectorized = solverSpecs.Vectorized;
            end
        end

        % % sys: ODESys
        % function setSolverEvents(sys)
        %     % function v = solverEvts(t,y,p)
        %     % 
        %     % end
        %     % function [stop,ye,p] = solverEvtCallbacks(te,ye,ie,p)
        %     %     stop = false;
        %     % 
        %     % end
        % 
        %     solverEvtsText = '@(t,y,p) [';
        %     solverEvtCallbacks_ye_Text = {};
        %     solverEvtCallbacks_p_Text = {};
        %     solverEvtCallbacksText = '@(te,ye,ie,p) deal(false,';
        %     evt_dir = [];
        %     if sys.IO_flowrates.L_o == 0
        %         L_o_evt_txt = '';
        %         L_o_callback_ye_txt = {};
        %         L_o_callback_p_txt = {};
        %     else
        %         % descending
        %         L_o_evt_txt = ['y(',char(string(sys.IO_flowrates.V_idx)),') - 0.01*',char(string(sys.PBR_max_vol))];
        %         L_o_callback_ye_txt = {};
        %         L_o_callback_p_txt = struct('',['(ie==1)*0+~(ie==1)*','p(',,')']);
        %         evt_dir(end+1) = "descending";
        %         % ascending
        %     end
        %     solverEvtCallbacks_ye_Text = [solverEvtCallbacks_ye_Text,{L_o_callback_ye_txt}];
        %     solverEvtCallbacks_p_Text = [solverEvtCallbacks_p_Text,{L_o_callback_p_txt}];
        %     if sys.IO_flowrates.L_i == 0
        %         L_i_evt_txt = '';
        %         L_i_callback_txt = '';
        %     else
        %         L_i_evt_txt = '';
        %         L_i_callback_txt = '';
        %         evt_dir(end+1) = "ascending";
        %     end
        %     solverEvtsText = [solverEvtsText,L_o_evt_txt];
        %     evt_response = repmat("callback",size(evt_dir));
        % end

        function generatePlots(sys)
            % running model for each plot
            sys.compileModel();

            % creating separate figures for each subplot
            sys.subplots = {};
            sys.TLs = {};
            for k=1:1:sys.subplot_ct
                sys.subplots{k} = figure(k);
                sys.TLs{k} = tiledlayout(sys.subplots{k},sys.subplot_row_ct,sys.subplot_col_ct);
            end
            test2 = sys.subplot_ct
            test3 = length(sys.plots)
            
            % ### FIXME: add feature to allow user to specify time
            % precision of model
            tEnd = sys.model_runtime;
            tSmooth = [0,tEnd];
            sysVar = sys.getModelVarNames("plot");
            for k=1:1:length(sys.plots)
                plot_obj = sys.plots{k};
                axes = plot_obj.axes;
                if plot_obj.getPlotProp("display") == true || plot_obj.getPlotProp("download") == true
                    if length(axes) == 3
                        if axes{1}.varIsIC && axes{2}.varIsIC
                            y0_span_x = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts);
                            y0_span_y = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts);
                            % maybe remove these and just use y0_mesh_X and
                            % y0_mesh_Y for the interpolation query points
                            % as well?
                            y0_span_qx = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts*5);
                            y0_span_qy = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts*5);
                            [y0_mesh_X,y0_mesh_Y] = meshgrid(y0_span_x,y0_span_y);
                            [y0_mesh_qX,y0_mesh_qY] = meshgrid(y0_span_qx,y0_span_qy);

                            res_1 = cell(length(y0_span_y),length(y0_span_x));

                            % creating array that maps y0 to variable name
                            y0_2_name = [];
                            y0_other_phase = [];
                            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                            ext_ct = 0;
                            for l=1:1:length(comps)
                                y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Liquid Phase)']); %#ok<AGROW>
                                if comps{l}.is_vol
                                    ext_ct = ext_ct + 1;
                                    % helperName = ['Gas Phase Partial Pressure of ',char(comps{l}.getName()),' in Equilibrium with Liquid Phase'];
                                    % [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                    % eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                    % % ### FIXME?: unable to access initial conditions for
                                    % % variables inputted from Simulink, will assume 0 for
                                    % % all
                                    % sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));

                                    y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Gas Phase)']); %#ok<AGROW>
                                    y0_other_phase(end+1) = l; %#ok<AGROW>
                                end
                                solvents = {};
                                if ~isempty(comps{k}.sorpFuncParams)
                                    for m=1:1:length(comps{k}.sorpFuncParams)
                                        if ~any(strcmp(solvents,comps{k}.sorpFuncParams{l}.solventName))
                                            ext_ct = ext_ct + 1;
                                            % helperName = [char(comps{k}.name),' Concentration in ',comps{k}.sorpFuncParams{l}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                            % [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            % eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            % sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));

                                            y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (',comps{k}.sorpFuncParams{l}.solventName,' Phase)']); %#ok<AGROW>
                                            y0_other_phase(end+1) = l; %#ok<AGROW>
                                            solvents{end+1} = comps{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                                        end
                                    end
                                end
                            end

                            y0_var_num_x = strcmp(y0_2_name,erase(axes{1}.getPlotProp("varNames"),"~(IC)"));
                            y0_var_num_y = strcmp(y0_2_name,erase(axes{2}.getPlotProp("varNames"),"~(IC)"));

                            multiple_x = length(find(y0_other_phase(y0_var_num_x) == y0_other_phase));
                            multiple_y = length(find(y0_other_phase(y0_var_num_y) == y0_other_phase));

                            for l=1:1:length(y0_span_y)
                                for m=1:1:length(y0_span_x)
                                    y0_1 = sys.y0;
                                    % ### FIXME: need to write logic to edit the
                                    % changes in initial conditions when there are
                                    % chemical solutes in different phases
                                    y0_1(y0_var_num_x) = y0_span_(m);
                                    if multiple_x
                                        first_x_idx = find(y0_other_phase(y0_var_num_x) == y0_other_phase,1);
                                        % setting comp
                                        comp = comps{y0_other_phase(first_x_idx)};
                                        if first_x_idx + 1 == y0_var_num_x
                                            helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            y0_1(first_x_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        elseif first_x_idx ~= y0_var_num_x
                                            for n=1:1:length(comp.sorpFuncParams)
                                                if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{n}.solventName)
                                                    solventName = comp.sorpFuncParams{n}.solventName;
                                                end
                                            end
                                            helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            y0_1(first_x_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        end

                                        % now update the other phases
                                        ext_ct = 0;
                                        if comp.is_vol
                                            ext_ct = ext_ct + 1;
                                            helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            % ### FIXME?: unable to access initial conditions for
                                            % variables inputted from Simulink, will assume 0 for
                                            % all
                                            y0_1(first_x_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        end
                                        solvents = {};
                                        if ~isempty(comp.sorpFuncParams)
                                            for n=1:1:length(comp.sorpFuncParams)
                                                if ~any(strcmp(solvents,comp.sorpFuncParams{n}.solventName))
                                                    ext_ct = ext_ct + 1;
                                                    helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{n}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                                    y0_1(first_x_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                                    solvents{end+1} = comp.sorpFuncParams{n}.solventName; %#ok<AGROW>
                                                end
                                            end
                                        end
                                    end
                                    y0_1(y0_var_num_y) = y0_span_y(l);
                                    if multiple_y
                                        first_y_idx = find(y0_other_phase(y0_var_num_y) == y0_other_phase,1);
                                        % setting comp
                                        comp = comps{y0_other_phase(first_y_idx)};
                                        if first_y_idx + 1 == y0_var_num_y
                                            helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            y0_1(first_y_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        elseif first_y_idx ~= y0_var_num_y
                                            for n=1:1:length(comp.sorpFuncParams)
                                                if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{n}.solventName)
                                                    solventName = comp.sorpFuncParams{n}.solventName;
                                                end
                                            end
                                            helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            y0_1(first_y_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        end

                                        % now update the other phases
                                        ext_ct = 0;
                                        if comp.is_vol
                                            ext_ct = ext_ct + 1;
                                            helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            % ### FIXME?: unable to access initial conditions for
                                            % variables inputted from Simulink, will assume 0 for
                                            % all
                                            y0_1(first_y_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                        end
                                        solvents = {};
                                        if ~isempty(comp.sorpFuncParams)
                                            for n=1:1:length(comp.sorpFuncParams)
                                                if ~any(strcmp(solvents,comp.sorpFuncParams{n}.solventName))
                                                    ext_ct = ext_ct + 1;
                                                    helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{n}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                                    y0_1(first_y_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                                    solvents{end+1} = comp.sorpFuncParams{n}.solventName; %#ok<AGROW>
                                                end
                                            end
                                        end
                                    end

                                    [tRes,yRes] = sys.runModel(tSmooth,y0_1);
                                    fRes = sys.calculateHelperVals(tRes,yRes);
                                    res_1{l,m} = [yRes,fRes,tRes];
                                end
                            end

                            res_2 = {};
                            evalt = sys.plots{k}.getPlotProp("evaltVal");
                            for l=1:1:length(axes{3}.getPlotProp("varNames"))
                                yVarIdx = strcmp(sysVar(:,1),axes{3}.varNames{l});
                                res_2{l} = zeros(size(res_1)); %#ok<AGROW>
                                for m=1:1:size(res_2,1)
                                    for n=1:1:size(res_2,2)
                                        % ### FIXME: add more methods for
                                        % interpolation?
                                        res_2{l}(m,n) = interp1(res_1{m,n}(:,end),res_1{m,n}(:,yVarIdx),evalt,'makima');
                                    end
                                end
                            end
                        elseif axes{1}.varIsIC || axes{2}.varIsIC
                            if axes{1}.varIsIC
                                y0_span = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts);
                                ICax = 1;
                                % ICax_obj = axes{1};
                                % nonICAx = 2;
                                % nonICAx_obj = axes{2};
                            else
                                y0_span = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts);
                                ICax = 2;
                                % ICax_obj = axes{2};
                                % nonICAx = 1;
                                nonICax_obj = axes{1};
                            end

                            res_1 = cell(length(y0_span),1);

                            % creating array that maps y0 to variable name
                            y0_2_name = [];
                            y0_other_phase = [];
                            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                            ext_ct = 0;
                            for l=1:1:length(comps)
                                y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Liquid Phase)']); %#ok<AGROW>
                                if comps{l}.is_vol
                                    ext_ct = ext_ct + 1;
                                    y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Gas Phase)']); %#ok<AGROW>
                                    y0_other_phase(end+1) = l; %#ok<AGROW>
                                end
                                solvents = {};
                                if ~isempty(comps{k}.sorpFuncParams)
                                    for m=1:1:length(comps{k}.sorpFuncParams)
                                        if ~any(strcmp(solvents,comps{k}.sorpFuncParams{l}.solventName))
                                            ext_ct = ext_ct + 1;
                                            y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (',comps{k}.sorpFuncParams{l}.solventName,' Phase)']); %#ok<AGROW>
                                            y0_other_phase(end+1) = l; %#ok<AGROW>
                                            solvents{end+1} = comps{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                                        end
                                    end
                                end
                            end

                            y0_var_num = strcmp(y0_2_name,erase(plot_obj.getAxProp(ICax,"varNames"),"~(IC)"));

                            multiple_num = length(find(y0_other_phase(y0_var_num) == y0_other_phase));

                            for l=1:1:length(y0_span)
                                y0_1 = sys.y0;
                                y0_1(y0_var_num) = y0_span(l);
                                if multiple_num
                                    first_idx = find(y0_other_phase(y0_var_num) == y0_other_phase,1);
                                    % setting comp
                                    comp = comps{y0_other_phase(first_idx)};
                                    if first_idx + 1 == y0_var_num
                                        helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    elseif first_idx ~= y0_var_num
                                        for m=1:1:length(comp.sorpFuncParams)
                                            if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{m}.solventName)
                                                solventName = comp.sorpFuncParams{m}.solventName;
                                            end
                                        end
                                        helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end

                                    % now update the other phases
                                    ext_ct = 0;
                                    if comp.is_vol
                                        ext_ct = ext_ct + 1;
                                        helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end
                                    solvents = {};
                                    if ~isempty(comp.sorpFuncParams)
                                        for m=1:1:length(comp.sorpFuncParams)
                                            if ~any(strcmp(solvents,comp.sorpFuncParams{m}.solventName))
                                                ext_ct = ext_ct + 1;
                                                helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{m}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                                [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                                eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                                y0_1(first_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                                solvents{end+1} = comp.sorpFuncParams{m}.solventName; %#ok<AGROW>
                                            end
                                        end
                                    end

                                    [tRes,yRes] = sys.runModel(tSmooth,y0_1);
                                    fRes = sys.calculateHelperVals(tRes,yRes);
                                    res_1{l,1} = [yRes,fRes,tRes];
                                end
                            end

                            res_2 = {};
                            nonICVarIdx = strcmp(sysVar(:,1),nonICax_obj.varNames);
                            for l=1:1:length(plot_obj.getAxProp(3,"varNames"))
                                zVarIdx = strcmp(sysVar(:,1),axes{3}.varNames{l});
                                res_2{l} = zeros(size(res_1)); %#ok<AGROW>
                                for m=1:1:size(res_2,1)
                                    % ### FIXME: add more methods for
                                    % interpolation?
                                    res_2{l}(m,1) = interp1(res_1{m,1}(:,nonICVarIdx),res_1{m,1}(:,zVarIdx), ...
                                        linspace(min(res_1{m,1}(:,nonICVarIdx)),max((res_1{m,1}(:,nonICVarIdx)),size(res_1{m,1},1))),'makima');
                                end
                            end
                        else
                            [tRes,yRes] = sys.runModel(tSmooth);
                            fRes = sys.calculateHelperVals(tRes,yRes);
                            res_1{1,1} = [yRes,fRes,tRes];

                            xVarIdx = strcmp(sysVar(:,1),axes{1}.varNames);
                            yVarIdx = strcmp(sysVar(:,1),axes{2}.varNames);
                            zVarIdx = zeros(size(axes{3}.varNames));
                            for l=1:1:length(axes{3}.varNames)
                                zVarIdx(l) = strcmp(sysVar(:,1),axes{2}.varNames{l});
                            end

                            y0_span_x = res_1{1,1}(:,xVarIdx);
                            y0_span_y = res_1{1,1}(:,yVarIdx);
                            y0_span_z = res_1{1,1}(:,zVarIdx);
                        end
                    else
                        [tRes,yRes] = sys.runModel(tSmooth,sys.y0);
                        fRes = sys.calculateHelperVals(tRes,yRes);
                        res_1{1,1} = [tRes,yRes,fRes];
                    end

                    % plot models on fig
                    % ### IMPROVEMENT: give user option to span multiple
                    % slots?
                    fig = figure(sys.plots{k}.subplotGroup);
                    group = sys.plots{k}.subplotGroup;
                    slot = sys.plots{k}.subplotSlot;
                    ax = nexttile(sys.TLs{group},slot);
                    hold(ax,"on");
                    if length(axes) == 2
                        xVarIdx = strcmp(sysVar(:,1),axes{1}.varNames);
                        yVarIdx = zeros(size(axes{2}.varNames));
                        for l=1:1:length(axes{2}.varNames)
                            yVarIdx(l) = find(strcmp(sysVar(:,1),axes{2}.varNames{l}));
                        end
                        plot(res_1{1,1}(:,xVarIdx),res_1{1,1}(:,yVarIdx));
                    elseif length(axes) == 3          
                        if axes{1}.varIsIC && axes{2}.varIsIC
                            for l=1:1:length(res_2)
                                % 2D interpolation
                                interpZGrid = interp2(y0_mesh_X, y0_mesh_Y, res_2{l}, y0_mesh_qX, y0_mesh_qY, 'makima');
                                % plotting surface
                                surf(qGridXVals,qGridYVals,interpZGrid);
                            end
                        elseif axes{1}.varIsIC || axes{2}.varIsIC
                            for l=1:1:length(res_2)
                                % 2D interpolation
                                interpZGrid = interp2(y0_mesh_X, y0_mesh_Y, res_2{l}, y0_mesh_qX, y0_mesh_qY, 'makima');
                                % plotting surface
                                surf(qGridXVals,qGridYVals,interpZGrid);
                            end
                        else
                            % plotting lines in 3D
                            plot3(y0_span_x,y0_span_y,y0_span_z);
                        end
                    end

                    title(plot_obj.title);
                    for l=1:1:length(plot_obj.axes)
                        if l == 1
                            xlabel(axes{l}.title);
                            if ~axes{l}.useDefR, xlim([axes{l}.loDispLim,axes{l}.upDispLim]); end
                        elseif l == 2
                            ylabel(axes{l}.title);
                            if ~axes{l}.useDefR, ylim([axes{l}.loDispLim,axes{l}.upDispLim]); end
                        elseif l == 3
                            zlabel(axes{l}.title);
                            if ~axes{l}.useDefR, zlim([axes{l}.loDispLim,axes{l}.upDispLim]); end
                        end
                    end
                    legend(axes{end}.varNames);
                    hold(ax,"off");

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
        function [tRes,yRes] = runModel(sys,tspan,y0)
            try
                % ### QUALITY: build in timeout error with events for ODE
                % solver
                sys.solver.ODEFun = @(t,y,p) sys.dydt(t,y,p,sys.f(:,3),zeros(size(sys.var_in_key)));
                sys.solver.InitialTime = 0;
                sys.solver.InitialValue = y0;
                sys.solver.Parameters = sys.param;
                S = solve(sys.solver,tspan(1),tspan(end));
                tRes = S.Time;
                yRes = S.Solution;
            catch err
                % ### ACTUALLYSTARTHERE: if there are errors in compiliation,
                % then throw them and alert the user; that way, it is up to the
                % user to modify their fermentation conditions to avoid
                % overflow or tank draining conditions
            end
        end

        % sys: ODESys class ref, tRes: num[], yRes: num[]
        function fRes = calculateHelperVals(sys,tRes,yRes)
            fRes = zeros(size(tRes,1),size(sys.f,1));
            for k=1:1:size(sys.f,1)
                for l=1:1:size(tRes,1)
                    fRes(l,k) = sys.f{k,3}(tRes(l),yRes(l,:),sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                end
            end
        end

        % sys: ODESys class ref,
        function regStats = compileRegression(sys)
            sys.dydt = "@(t,y,p,f,v,r) [";
            sys.param = [];
            sys.reg_param_ct = 0;
            sys.f = {};
            for k=1:1:size(sys.regParamList,1), sys.regParamList{k,6} = string(k); end

            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
            sys.var_out_key = cell(length(comps),2);
            sys.S_nums = {};
            ext_ct = 0;
            for k=1:1:length(comps)
                sys.var_out_key{k+ext_ct,1} = comps{k}.getName();
                sys.var_out_key{k+ext_ct,2} = comps{k}.getSym();
                if comps{k}.is_vol  
                    ext_ct = ext_ct + 1;
                    sys.var_out_key{k+ext_ct,2} = comps{k}.getName();
                    sys.var_out_key{k+ext_ct,2} = comps{k}.getBulkGasSym();
                end
                for l=1:1:length(comps{k}.sorpFuncParams)
                    ext_ct = ext_ct + 1;
                    sys.var_out_key{k+ext_ct,2} = comps{k}.getName();
                    sys.var_out_key{k+ext_ct,2} = comps{k}.sorpFuncParams{l}.solventSym;
                end
                if strcmp(comps{k}.getType(),'Suspended Solid Sorbent')
                    sys.S_nums{end+1} = comps{k}.getNum();
                end
            end

            in_syms = sys.getModelVarNames("updateModel");
            in_syms = in_syms(:,2);
            sys.var_in_key = {};
            num_spec = length(sys.getSpecies('comp'));
            num_spec_chem = length([sys.getSpecies('comp'),sys.getChemicals('comp')]);

            sysVar = sys.getModelVarNames("plot");
            subfuncs = [sys.environs.(sys.activeEnv).subfuncs,sys.helperFuncs];
            for k=1:1:length(comps)
                % add logic to compile component gov funcs in order that
                % allows for growth-associated product funcs to be a
                % function of biomass funcs, and substrate funcs to be a
                % function of biomass and product funcs
                % ### FIXME: requires special logic/restrictions in the
                % governing function portion
                [govFunc, params, sys.reg_param_ct, regPUpdate] = comps{k}.compileGovFunc(length(sys.param),sys.reg_param_ct,sys.regParamList(:,[1,6]),true);
                for l=1:1:size(regPUpdate,2)
                    if regPUpdate{1,l}, sys.regParamList{regPUpdate{2,l},6} = regPUpdate{3,l}; end
                end
                for l=1:1:length(subfuncs)
                    if ~strcmp(subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,subfuncs{l}.getSubFuncSym(),"f{"+l+"}(t,y,p,f,v,r)");
                    end
                end
                for l=1:1:length(in_syms)
                    if ~isempty(regexp(govFunc,in_syms{l},'once'))
                        sys.var_in_key{end+1} = in_syms{l};
                        govFunc = replace(govFunc,in_syms{l},['v(',char(string(length(sys.var_in_key))),')']);
                    end
                end

                ext_ct = 0;
                for l=1:1:length(comps)
                    if l <= num_spec
                        govFunc = regexprep(govFunc,"X_"+string(comps{l}.getNum()),"y("+l+")");
                    elseif l > num_spec && l <= num_spec_chem
                        govFunc = replace(govFunc,['C_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_C_',char(string(comps{l}.getNum()))]);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['C_',char(string(comps{l}.getNum())),'_']);
                        if comps{l}.is_vol
                            ext_ct = ext_ct + 1;
                            govFunc = regexprep(govFunc,"P_"+(comps{l}.getNum()),"y("+string(l+ext_ct)+")");
                        end
                        phases = comps{l}.getSorpPhases();
                        for m=1:1:length(phases)
                            ext_ct = ext_ct + 1;
                            govFunc = regexprep(govFunc,phases{m},"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_',phases{m}]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],[phases{m},'_']);
                        end
                        if strcmp(comps{l}.getType(),'Suspended Solid Sorbent')
                            govFunc = replace(govFunc,['V_S_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_V_S_',char(string(comps{l}.getNum()))]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['V_S_',char(string(comps{l}.getNum())),'_']);
                        end
                    elseif l == num_spec_chem + 1
                        govFunc = regexprep(govFunc,'T',['y(',char(string(l+ext_ct)),')']);
                    elseif l == num_spec_chem + 2
                        govFunc = regexprep(govFunc,'P',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'P_');
                    elseif l == num_spec_chem + 3
                        govFunc = regexprep(govFunc,'V',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'V_');
                    elseif l == num_spec_chem + 4
                        govFunc = replace(govFunc,'H3O',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_H3O');
                    elseif l == num_spec_chem + 5
                        govFunc = replace(govFunc,'OH',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_OH');
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
            for k=1:1:length(subfuncs)
                % ### FIXME: add feature to define helpers with another
                % helper (this is very complicated, for now just don't let
                % helpers be defined by helpers)
                ext_ct = 0;
                [govFunc,params,regPUpdate] = subfuncs{k}.compileSubFunc(length(sys.param),sys.reg_param_ct,sys.regParamList(:,[1,6]),true);
                test19 = size(regPUpdate,2)
                if size(regPUpdate,1) > 1
                    for l=1:1:size(regPUpdate,2)
                        if regPUpdate{1,l}, sys.regParamList{regPUpdate{2,l},6} = regPUpdate{3,l}; end
                    end
                end
                for l=1:1:length(comps)
                    if l <= num_spec
                        govFunc = replace(govFunc,"X_"+string(comps{l}.getNum()),"y("+l+")");
                    elseif l > num_spec && l <= num_spec_chem
                        govFunc = replace(govFunc,['C_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_C_',char(string(comps{l}.getNum()))]);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['C_',char(string(comps{l}.getNum())),'_']);
                        if comps{l}.is_vol
                            ext_ct = ext_ct + 1;
                            govFunc = replace(govFunc,"P_"+(comps{l}.getNum()),"y("+string(l+ext_ct)+")");
                        end
                        phases = comps{l}.getSorpPhases();
                        for m=1:1:length(phases)
                            ext_ct = ext_ct + 1;
                            govFunc = replace(govFunc,phases{m},"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_',phases{m}]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],[phases{m},'_']);
                        end
                        if strcmp(comps{l}.getType(),'Suspended Solid Sorbent')
                            govFunc = replace(govFunc,['V_S_',char(string(comps{l}.getNum()))],"y("+(l+ext_ct)+")");
                            govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],['_V_S_',char(string(comps{l}.getNum()))]);
                            govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],['V_S_',char(string(comps{l}.getNum())),'_']);
                        end
                    elseif l == num_spec_chem + 1
                        govFunc = regexprep(govFunc,'T',['y(',char(string(l+ext_ct)),')']);
                    elseif l == num_spec_chem + 2
                        govFunc = regexprep(govFunc,'P',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'P_');
                    elseif l == num_spec_chem + 3
                        govFunc = regexprep(govFunc,'V',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['y\(',char(string(l+ext_ct)),'\)_'],'V_');
                    elseif l == num_spec_chem + 4
                        govFunc = replace(govFunc,'H3O',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_H3O');
                    elseif l == num_spec_chem + 5
                        govFunc = replace(govFunc,'OH',['y(',char(string(l+ext_ct)),')']);
                        govFunc = regexprep(govFunc,['_y\(',char(string(l+ext_ct)),'\)'],'_OH');
                    end
                end
                for l=1:1:length(subfuncs)
                    if ~strcmp(subfuncs{l}.getSubFuncSym(),"t")
                        govFunc = regexprep(govFunc,subfuncs{l}.getSubFuncSym(),"f{"+l+"}(t,y,p,f,v,r)");
                    end
                end
                % making sure the in_syms are tracked accurately from
                % component govFuncs
                for l=1:1:length(in_syms)
                    if ~isempty(regexp(govFunc,in_syms{l},'once'))
                        if any(strcmp(sys.var_in_key,in_syms{l}))
                            govFunc = replace(govFunc,in_syms{l},['v(',char(string(find([sys.var_in_key{:}] == in_syms{l}))),')']);
                        else
                            sys.var_in_key{end+1} = in_syms{l};
                            govFunc = replace(govFunc,in_syms{l},['v(',char(string(length(sys.var_in_key))),')']);
                        end
                    end
                end

                funcArgText = "@(t,y,p,f,v,r)";
                sys.f{end+1,1} = subfuncs{k}.getSubFuncName();
                sys.f{end,2} = subfuncs{k}.getSubFuncSym();
                sys.f{end,3} = str2func(funcArgText+govFunc);
                sys.param = [sys.param,params];
            end

            % ### STARTHERE: calculate values for v() array
            comps_v = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            for k=1:1:length(sys.var_in_key)
                switch sys.var_in_key{k}
                    case 'dP_o'


                    case 'L_i'


                    case 'L_o'


                    case 'D_i'

                    otherwise
                        if contains(sys.var_in_key{k},'L_s_i_')

                        elseif contains(sys.var_in_key{k},'L_s_o_')

                        elseif contains(sys.var_in_key{k},'D_i_')

                        elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'C_') && ~contains(sys.var_in_key{k},'D_') && ~contains(sys.var_in_key{k},'L_s_')

                        elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'Y_')

                        elseif contains(sys.var_in_key{k},'_i') && contains(sys.var_in_key{k},'S_') && ~contains(sys.var_in_key{k},'D_') && ~contains(sys.var_in_key{k},'L_s_')

                        end
                end
            end

            % converting to function_handle
            sys.dydt = str2func(sys.dydt');

            test = sys.dydt
            test5 = sys.f{1,3}
            test2 = sys.f{7,3}

            sys.importedDataIdx = [];
            for k=1:1:length(sys.matchedVarsList)
                matchIdx = strcmp(cellstr(sysVar(:,1)),sys.matchedVarsList{k}.sysVarName);
                varSyms = cellstr(sysVar(:,1));
                if any(matchIdx) && ~strcmp("Model Runtime",varSyms{matchIdx})
                    sys.importedDataIdx(end+1) = find(matchIdx);
                end
            end
            sys.IVs = zeros((size(sys.importedData,2)-1).*size(sys.importedData,1),1);
            sys.DVs = zeros((size(sys.importedData,2)-1).*size(sys.importedData,1),1);
            sys.mat_IVs = [sys.importedData{:,1}];
            sys.mat_DVs = zeros(size(sys.importedData,1),size(sys.importedData,2)-1);
            for k=2:1:size(sys.importedData,2)
                sys.IVs((1:1:size(sys.importedData,1))+((k-2).*size(sys.importedData,1)),1) = [sys.importedData{:,1}];
                sys.mat_DVs(:,k-1) = sys.importedData{:,k};
                % if strcmp(comps{sys.importedDataIdx(k-1)-1}.getInitConcUnit(),"mg/L")
                %     sys.DVs((1:1:size(sys.importedData,1))+((k-2).*size(sys.importedData,1)),1) = sys.importedData{:,k}./1000;
                % else
                %     sys.DVs((1:1:size(sys.importedData,1))+((k-2).*size(sys.importedData,1)),1) = sys.importedData{:,k};
                % end
                sys.DVs((1:1:size(sys.importedData,1))+((k-2).*size(sys.importedData,1)),1) = unit_standardization(sys.importedData{:,k},comps{sys.importedDataIdx(k-1)-1}.getInitConcUnit());
            end
            % ### FIXME: upgrade later to allow automatic testing of
            % various starting guesses for each parameter
            beta0 = sys.regSpecs.paramIGs;

            sys.setInitCond();

            % set proper reg specs based on solver
            switch sys.regSpecs.solver
                case "nlinfit"
                    nlinfitSpecs = statset;
                    nlinfitSpecs.Display = 'iter';
                    nlinfitSpecs.MaxIter = sys.regSpecs.nlinfitSpecs.MaxIter;
                    nlinfitSpecs.TolFun = sys.regSpecs.nlinfitSpecs.TolFun;
                    nlinfitSpecs.TolX = sys.regSpecs.nlinfitSpecs.TolX;
                    nlinfitSpecs.DerivStep = sys.regSpecs.nlinfitSpecs.DerivStep;
                    nlinfitSpecs.FunValCheck = sys.regSpecs.nlinfitSpecs.FunValCheck;
                    nlinfitSpecs.RobustWgtFun = sys.regSpecs.nlinfitSpecs.RobustWgtFun;
                    if sys.regSpecs.nlinfitSpecs.Tune >= 0, nlinfitSpecs.Tune = sys.regSpecs.nlinfitSpecs.Tune; end

                    [sys.reg_analytics.beta,sys.reg_analytics.R,sys.reg_analytics.J, ...
                        sys.reg_analytics.CovB,sys.reg_analytics.MSE,sys.reg_analytics.ErrorModelInfo] = ...
                        nlinfit(sys.IVs,sys.DVs,@(reg_param,t) sys.nlinfitHandler(reg_param,t,false),beta0,nlinfitSpecs);

                case "lsqcurvefit"
                    lsqcurvefitSpecs = optimset;
                    lsqcurvefitSpecs.Display = 'iter';
                    lsqcurvefitSpecs.Algorithm = sys.regSpecs.lsqcurvefitSpecs.Algorithm;
                    lsqcurvefitSpecs.PlotFcns = sys.regSpecs.lsqcurvefitSpecs.PlotFcns;
                    lsqcurvefitSpecs.TolFun = sys.regSpecs.lsqcurvefitSpecs.TolFun;
                    lsqcurvefitSpecs.TolX = sys.regSpecs.lsqcurvefitSpecs.TolX;
                    lsqcurvefitSpecs.MaxIter = sys.regSpecs.lsqcurvefitSpecs.MaxIter;
                    lsqcurvefitSpecs.MaxFunEvals = sys.regSpecs.lsqcurvefitSpecs.MaxFunEvals;
                    lsqcurvefitSpecs.FunValCheck = sys.regSpecs.lsqcurvefitSpecs.FunValCheck;
                    [sys.reg_analytics.beta,~,sys.reg_analytics.R,~,~,~,sys.reg_analytics.J] = ...
                        lsqcurvefit(@(reg_param,t) sys.nlinfitHandler(reg_param,t,false),beta0,sys.IVs,sys.DVs, ...
                        sys.regSpecs.lsqcurvefitSpecs.paramBnds(:,1),sys.regSpecs.lsqcurvefitSpecs.paramBnds(:,2), ...
                        [],[],[],[],[],lsqcurvefitSpecs);
                    sys.reg_analytics.MSE = mean(sys.reg_analytics.R.^2,"all");

                case "fminsearch"
                    fminsearchSpecs = optimset;
                    fminsearchSpecs.Display = 'iter';
                    fminsearchSpecs.FunValCheck = sys.regSpecs.fminsearchSpecs.FunValCheck;
                    fminsearchSpecs.MaxFunEvals = sys.regSpecs.fminsearchSpecs.MaxFunEvals;
                    fminsearchSpecs.MaxIter = sys.regSpecs.fminsearchSpecs.MaxIter;
                    fminsearchSpecs.PlotFcns = sys.regSpecs.fminsearchSpecs.PlotFcns;
                    fminsearchSpecs.TolFun = sys.regSpecs.fminsearchSpecs.TolFun;
                    fminsearchSpecs.TolX = sys.regSpecs.fminsearchSpecs.TolX;
                    sys.reg_analytics.beta = fminsearch(@(reg_param) sys.fminsearchHandler(sys.IVs,reg_param),beta0,fminsearchSpecs);
                    sys.reg_analytics.R = (sys.DVs-sys.runRegModel(sys.IVs,sys.y0,sys.reg_analytics.beta));
                    % sys.reg_analytics.MSE = mean((DVs-sys.runRegModel(IVs,sys.y0,sys.reg_analytics.beta)).^2,"all");
                    sys.reg_analytics.MSE = mean(sys.reg_analytics.R.^2,"all");

                case "fmincon"
                    fminconSpecs = optimset;
                    fminconSpecs.Algorithm = sys.regSpecs.fminconSpecs.Algorithm;
                    fminconSpecs.PlotFcns = sys.regSpecs.fminconSpecs.PlotFcns;
                    fminconSpecs.TolFun = sys.regSpecs.fminconSpecs.TolFun;
                    fminconSpecs.TolX = sys.regSpecs.fminconSpecs.TolX;
                    fminconSpecs.MaxIter = sys.regSpecs.fminconSpecs.MaxIter;
                    fminconSpecs.MaxFunEvals = sys.regSpecs.fminconSpecs.MaxFunEvals;
                    fminconSpecs.FunValCheck = sys.regSpecs.fminconSpecs.FunValCheck;
                    sys.reg_analytics.beta = lsqcurvefit(@(reg_param) sys.fminsearchHandler(sys.IVs,reg_param),beta0, ...
                            [],[],[],[],sys.regSpecs.fminconSpecs.paramBnds(:,1),sys.regSpecs.fminconSpecs.paramBnds(:,2), ...
                            [],fminconSpecs);
                    sys.reg_analytics.R = (sys.DVs-sys.runRegModel(sys.IVs,sys.y0,sys.reg_analytics.beta));
                    sys.reg_analytics.MSE = mean((sys.DVs-sys.runRegModel(sys.IVs,sys.y0,sys.reg_analytics.beta)).^2,"all");

            end

            [tRes,yRes] = sys.runRegModel(sys.IVs,sys.y0,sys.reg_analytics.beta);
            test76 = sys.reg_analytics.beta

            sys.regData = [tRes,yRes];
            for k=1:1:length(sys.importedDataIdx)
                sys.regData(:,sys.importedDataIdx(k)) = unit_standardization_inv(sys.regData(:,sys.importedDataIdx(k)),comps{sys.importedDataIdx(k)-1}.getInitConcUnit());
                % if strcmp(comps{sys.importedDataIdx(k)-1}.getInitConcUnit(),"mg/L")
                %     sys.regData(:,sys.importedDataIdx(k)) = sys.regData(:,sys.importedDataIdx(k)).*1000;
                % end
            end
            regStats = struct('importedDataIdx',sys.importedDataIdx, ...
                'beta',sys.reg_analytics.beta,'R',sys.reg_analytics.R, ...
                'J',sys.reg_analytics.J,'CovB',sys.reg_analytics.CovB,'MSE',sys.reg_analytics.MSE, ...
                'ErrorModelInfo',sys.reg_analytics.ErrorModelInfo, ...
                'IVs',sys.mat_IVs,'DVs',sys.mat_DVs,'tRes',sys.regData(:,1),'yRes',sys.regData(:,2:end));
        end

        % sys: ODESys class ref, param: number[], t: number[]
        function yRes = nlinfitHandler(sys,reg_param,t,varargin)
            % test = reg_param
            if varargin{1}, tspan = [0,t]; else, tspan = unique(t); end
            try
                sys.solver.ODEFun = @(t,y,p) sys.dydt(t,y,p,sys.f(:,3),zeros(size(sys.var_in_key)),reg_param);
                sys.solver.InitialTime = 0;
                sys.solver.InitialValue = sys.y0;
                sys.solver.Parameters = sys.param;
                S = solve(sys.solver,tspan(1),tspan(end));
                yRes_all = S.Solution;
                % [~,yRes_all] = ode15s(@(t,y) sys.dydt(t,y,sys.param,sys.f(:,3),zeros(size(sys.var_in_key,2)),reg_param),tspan,sys.y0,opts);
            catch err

            end
            if varargin{1}
                yRes = yRes_all(end,sys.importedDataIdx-1)';
            else
                yRes = zeros((size(sys.importedData,2)-1).*size(sys.importedData,1),1);
                for k=1:1:size(sys.importedData,2)-1
                    yRes((1:1:size(sys.importedData,1))+((k-1).*size(sys.importedData,1)),1) = yRes_all(:,sys.importedDataIdx(k)-1);
                end
                % test2 = yRes
                % test3 = reg_param
            end
            % test2 = yRes
        end

        % sys: ODESys class ref, x: number[]
        function SSE = fminsearchHandler(sys,IVs,reg_param)
            tspan = unique(IVs);
            try
                sys.solver.ODEFun = @(t,y,p) sys.dydt(t,y,p,sys.f(:,3),zeros(size(sys.var_in_key)),reg_param);
                sys.solver.InitialTime = 0;
                sys.solver.InitialValue = sys.y0;
                sys.solver.Parameters = sys.param;
                S = solve(sys.solver,tspan(1),tspan(end));
                yRes_all = S.Solution;
            catch err
                
            end
            % [~,yRes_all] = ode15s(@(t,y) sys.dydt(t,y,sys.param,sys.f(:,3),zeros(size(sys.var_in_key,2)),reg_param),tspan,sys.y0,opts);
            yRes = zeros(size(table2array(sys.importedData(:,2:end))));
            for k=1:1:size(sys.importedData,2)-1
                yRes(:,k) = yRes_all(:,sys.importedDataIdx(k));
            end
            % test3 = (yRes)
            % test4 = (table2array(sys.importedData(:,2:end)))
            SSE = sum((table2array(sys.importedData(:,2:end))-yRes).^2,"all");
        end

        % sys: ODESys class ref, t: number[]
        function [tRes,yRes] = runRegModel(sys,t,y0,reg_param)
            % iterating over BatchFunction with ode15s
            tspan = unique(t);
            try
                sys.solver.ODEFun = @(t,y,p) sys.dydt(t,y,p,sys.f(:,3),zeros(size(sys.var_in_key)),reg_param);
                sys.solver.InitialTime = 0;
                sys.solver.InitialValue = y0;
                sys.solver.Parameters = sys.param;
                S = solve(sys.solver,tspan(1),tspan(end));
                tRes = S.Time;
                yRes = S.Solution;
            catch err

            end
            % [tRes,yRes] = ode15s(@(t,y) sys.dydt(t,y,sys.param,sys.f(:,3),zeros(size(sys.var_in_key,2)),reg_param),tspan,y0,opts);
        end

        % sys: ODESys class ref
        function reg_analytics = getRegAnalytics(sys)
            reg_analytics = sys.reg_analytics;
        end

        % sys: ODESys class ref
        function NRMSE = calcNRMSE(sys,normType)
            % ### FIXME: check that calculation is correct
            RMSE = zeros(size(sys.importedData,2)-1,1);
            NRMSE = cell(size(sys.importedData,2)-1,3);
            skip = 0;
            for k=1:1:length(RMSE)
                RMSE(k) = sqrt(mean(sys.reg_analytics.R((1:1:size(sys.importedData,1))+((k-1).*size(sys.importedData,1)),1).^2));
                switch normType
                    case "Range"
                        normFactor = range(sys.regData(:,k+1),'all');
                    case "Mean"
                        normFactor = mean(sys.regData(:,k+1),'all');
                    case "Inter-Quartile Range"
                        normFactor = iqr(sys.regData(:,k+1),'all');
                    case "Median"
                        normFactor = median(sys.regData(:,k+1),'all');
                end
                % test = sys.matchedVarsList{k}.sysVarName
                if strcmp(sys.matchedVarsList{k}.sysVarName,'Model Runtime'), skip = 1; end
                NRMSE{k,1} = char(sys.matchedVarsList{k+skip}.sysVarName);
                NRMSE{k,2} = RMSE(k)./normFactor;
                NRMSE{k,3} = 1 - sum((sys.regData(:,sys.importedDataIdx(k))-sys.importedData{:,k+1}).^2)./sum((sys.importedData{:,k+1}-mean(sys.importedData{:,k+1})).^2);
            end
        end

        % sys: ODESys class ref, calcType: string, alpha: number
        function parCI = calcParCI(sys,calcType,alpha)
            % ### FIXME: need to ensure that the R, J, and CovB are
            % separated by variable
            % ### FIXME: Check that the calculation here is correct
            parCI = cell(size(sys.regParamList,1),3);
            switch calcType
                case "Covariance"
                    beta = sys.reg_analytics.beta;
                    R = sys.reg_analytics.R;
                    CovB = sys.reg_analytics.CovB;
                    for k=1:1:length(beta)
                        res = nlparci(beta,R,'covariance',CovB,'alpha',alpha);
                        parCI{k,2} = res(k,1);
                        parCI{k,3} = res(k,2);
                    end
                case "Jacobian"
                    beta = sys.reg_analytics.beta;
                    R = sys.reg_analytics.R;
                    J = sys.reg_analytics.J;
                    for k=1:1:length(beta)
                        res = nlparci(beta,R,'jacobian',J,'alpha',alpha);
                        parCI{k,2} = res(k,1);
                        parCI{k,3} = res(k,2);
                    end
            end
            for k=1:1:size(parCI,1), parCI{k,1} = sys.regParamList{k,1}; end
        end

        % sys: ODESys class ref, tval: number, calcType: string, alpha:
        % number, predIntType: string
        function predCI = calcPredCI(sys,tval,calcType,alpha,simBnds,predIntType)
            predCI = cell(length(sys.matchedVarsList)-1,4);
            if predIntType == "Curve", predIntType = 'curve'; else, predIntType = 'observation'; end
            if simBnds == "On", simBnds = 'on'; else, simBnds = 'off'; end

            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
            switch calcType
                case "Covariance"
                    % res = nlpredci(@(reg_param,t) sys.nlinfitHandler(reg_param,t,true),[tval,tval.*1.001],sys.reg_analytics.beta, ...
                    %     sys.reg_analytics.R,'covariance',sys.reg_analytics.CovB,'alpha',alpha, ...
                    %     'ErrorModelInfo',sys.reg_analytics.ErrorModelInfo,'MSE',sys.reg_analytics.MSE, ...
                    %     'PredOpt',predIntType,'SimOpt',simBnds);
                    [res,delta] = nlpredci(@(reg_param,t) sys.nlinfitHandler(reg_param,t,true),[tval,tval.*1.1],sys.reg_analytics.beta, ...
                        sys.reg_analytics.R,'covariance',sys.reg_analytics.CovB,'alpha',alpha,'MSE',sys.reg_analytics.MSE, ...
                        'PredOpt',predIntType,'SimOpt',simBnds);
                case "Jacobian"
                    [res,delta] = nlpredci(@(reg_param,t) sys.nlinfitHandler(reg_param,t,true),[tval,tval.*1.1],sys.reg_analytics.beta, ...
                        sys.reg_analytics.R,'jacobian',sys.reg_analytics.J,'alpha',alpha, ...
                        'MSE',sys.reg_analytics.MSE, ...
                        'PredOpt',predIntType,'SimOpt',simBnds);
            end

            for k=1:1:size(predCI,1)
                if strcmp(sys.matchedVarsList{k}.sysVarName,'Model Runtime'), skip = 1; end
                predCI{k,1} = char(sys.matchedVarsList{k+skip}.sysVarName);
                predCI{k,2} = res(k)-delta(k);
                predCI{k,3} = res(k)+delta(k);
                predCI{k,4} = comps{sys.importedDataIdx(k)-1}.getInitConcUnit();
                % predCI{k,4} = char(sys.matchedVarsList{k+skip}.sysVarUnit);
            end
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
                    sys.environs.(sys.activeEnv).subfuncs{k}.updateParams(paramName,paramSym,paramVal,paramUnit,true,sys.getDefaultParamVals());
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
            [group,slot] = sys.getNextSubplotSlot();
            sys.subplot_ct = sys.subplot_ct + 1;
            valid_name_found = false;
            plot_names = string(zeros(size(sys.plots)));
            for k=1:1:length(sys.plots), plot_names(k) = sys.plots{k}.getPlotProp("title"); end
            plot_name_num = length(sys.plots);
            while ~valid_name_found
                plot_name_num = plot_name_num+1;
                if ~any(strcmp(plot_names,"Plot "+plot_name_num))
                    valid_name_found = true;
                end
            end
            sys.plots{end+1} = Plot("Plot "+(plot_name_num),sys.getModelVarNames("plot"),sys.getModelVarNames("plot"),group,slot);
            plot = sys.plots{end}.getAllPlotProps();
            axes = sys.plots{end}.getAllAxProps();
        end

        % sys: ODESys class ref
        function [plot,axes] = removePlot(sys,plotName,lastItem)
            sys.subplot_ct = sys.subplot_ct - 1;
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

        % sys: ODESys class ref
        function clearPlotDispDownload(sys)
            for k=1:1:length(sys.plots)
                sys.plots{k}.updatePlot("display",false);
                sys.plots{k}.updatePlot("download",false);
            end
        end

        % sys: ODESys class ref, plotName: string, group: number, slot: number
        function updatePlotGroupSlot(sys,plotName,group,slot)
            plot = sys.getPlotByName(plotName);
            plot.subplotGroup = group;
            plot.subplotSlot = slot;
        end

        % sys: ODESys class ref
        function [group,slot] = getNextSubplotSlot(sys)
            filled_slots = cell(sys.subplot_ct+1,1);
            for k=1:1:length(filled_slots)
                filled_slots{k} = [];
            end
            if isempty(sys.plots)
                group = 1;
                slot = 1;
            else
                for k=1:1:length(sys.plots)
                    [group,slot] = sys.plots{k}.getSubplotGroupAndSlot();
                    filled_slots{group} = [filled_slots{group},slot];
                end
                for k=1:1:length(filled_slots), filled_slots{k} = sort(filled_slots{k}); end
                for k=1:1:length(filled_slots)
                    if length(filled_slots{k}) < (sys.subplot_row_ct*sys.subplot_col_ct)
                        group = k;
                        non_intersects = setxor(linspace(1,sys.subplot_row_ct*sys.subplot_col_ct,sys.subplot_row_ct*sys.subplot_col_ct),filled_slots{k});
                        slot = non_intersects(1);
                    end
                end
    
                % if all slots are full, throw an error
            end
        end

        % sys: ODESys class ref, group: number
        function updateSubplotPropsinPlot(sys,group)
            for k=1:1:length(sys.plots)
                if sys.plots{k}.subplotGroup == group
                    sys.plots{k}.subplotGroup = 0;
                    sys.plots{k}.subplotSlot = 0;
                end
            end
        end

        % sys: ODESys class ref, subplotNum: number
        function [row,col] = getSubplotRowsAndCols(sys,subplotNum)
            row = sys.subplot_row_ct(subplotNum);
            col = sys.subplot_col_ct(subplotNum);
        end

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
        function axes = updateAxVars(sys,plotName,axesDir,varName,addVar)
            plot = sys.getPlotByName(plotName);
            if contains(varName,"~")
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

            axes = plot.getAllAxProps();

            % varText = "Variable(s): $";
            % varNameStrings = convertCharsToStrings(plot.getAxProp(axesDir,"varNames"));
            % for k=1:1:length(varNameStrings)
            %     if k == length(varNameStrings)
            %         varText = varText + varNameStrings(k) + "$";
            %         return;
            %     end
            %     varText = varText + varNameStrings(k) + ", ";
            % end
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
            % ### FIXME: add logic to determine if varIsIC
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
        % varName: string, varUnits: string, evaltVal: string | number, loEvalLim: string | number, 
        % upEvalLim: string | number, nbEvalPts: number
        function setVarICEvalData(sys,plotName,axisName,varName,varUnits,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            plot = sys.getPlotByName(plotName);
            plot.setVarICEvalData(axisName,varName,varUnits,evaltVal,loEvalLim,upEvalLim,nbEvalPts);
        end

        % sys: ODESys class ref, filePath: string, colNames: {}
        function setImportedData(sys,filePath,importedData,colNames)
            sys.importedDataPath = filePath;
            sys.importedData = importedData;
            sys.importedDataColNames = colNames;
        end

        % sys: ODESys class ref
        function [importedData,path,colNames] = getUserImportedData(sys)
            importedData = sys.importedData;
            path = sys.importedDataPath;
            colNames = sys.importedDataColNames;
        end

        % sys: ODESys class ref, compName: string, paramSym: string, updateType: string
        function params = updateRegParamList(sys,compName,paramSym,updateType)
            try
                % ### FIXME: add functionality to include capability to update
                % for envs and helpers
                % comp = sys.getCompByName(compName);
                paramList = sys.getGrthParamsByFuncName(compName);
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
                        new_idx = 0;
                        for k=1:1:size(sys.regParamList,1)
                            sorted_comp_arr = sort({sys.regParamList{k,1},regParam{3}});
                            if strcmp(sorted_comp_arr{1},regParam{3})
                                new_idx = k;
                                break;
                            end
                        end
    
                        if ~isempty(sys.regParamList) && new_idx ~= 0
                            sys.regParamList(new_idx+1:size(sys.regParamList,1)+1,:) = sys.regParamList(new_idx:size(sys.regParamList,1),:);
                            sys.regSpecs.paramIGs(new_idx+1:size(sys.regSpecs.paramIGs,1)+1,1) = sys.regSpecs.paramIGs(new_idx:size(sys.regSpecs.paramIGs,1),1);
                        end
    
                        if new_idx == 0
                            new_idx = size(sys.regParamList,1)+1;
                        end

                        val = unit_standardization(regParam{5},regParam{6});
    
                        sys.regParamList{new_idx,1} = regParam{3};
                        sys.regParamList{new_idx,2} = compName;
                        sys.regParamList{new_idx,3} = regParam{2};
                        sys.regParamList{new_idx,4} = val;
                        sys.regParamList{new_idx,5} = '~';
                        sys.regParamList{new_idx,6} = "";
                        sys.regParamList{new_idx,7} = regParam{6};
    
                        % sys.regSpecs.DerivStep(end+1,1) = eps^(1/3);
                        sys.regSpecs.paramIGs(new_idx,1) = regParam{5};
                    end
                elseif updateType == "Remove"
                    for k=1:1:size(sys.regParamList,1)
                        if strcmp(sys.regParamList{k,1},paramSym)
                            sys.regParamList(k,:) = [];
                            % sys.regSpecs.DerivStep(k,1) = [];
                            sys.regSpecs.paramIGs(k,:) = [];
                            break;
                        end
                    end
                end
                % ### FIXME: for each of the regression solvers, need to update the
                % number of parameter specific specifications provided
                for k=1:1:size(sys.regParamList,1), sys.regParamList{k,6} = ""; end
                params = sys.regParamList(:,1:5);
            catch err
                disp("ODESys.m, line 2702: " + err)
            end
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
            % comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
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
                    % ### FIXME: need get the unit from the comp
                    % for k=1:1:length(comps)
                    % 
                    % end
                    sysVarUnit = 'g/L';
                    newPair = struct('sysVarName',sysVar,'sysVarNum',sysVarNum,'importVarName',importVar,'importVarNum',importVarNum,'sysVarUnit',sysVarUnit);
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
            plotRegData(:,2) = sys.regData(:,find(DVIdx)+1);
                        
            plotImportData = [];
            plotImportData(:,1) = [sys.importedData{:,1}];
            plotImportData(:,2) = [sys.importedData{:,find(sys.importedDataIdx-1 == (find(DVIdx)))+1}];
        end

        % sys: ODESys class ref
        function regSpecs = getRegSpecs(sys)
            regSpecs = sys.regSpecs;
        end

        % sys: ODEsys class ref, regSpecs: struct
        function setRegSpecs(sys,regSpecs)
            sys.regSpecs = regSpecs;
        end

        % sys: ODESys class ref, compName: string, phaseA: string, phaseB:
        % string
        function currRawGovFunc = getMTGovFunc(sys,compName,phaseA,phaseB)
            comp = sys.getCompByName([compName,' (Liquid Phase)']);
            currRawGovFunc = comp.getMTGovFunc(phaseA,phaseB);
        end

        % sys: ODESys class ref, funcName: string
        function MTFunc = getPresetMTFuncByName(sys,funcName,compName,phaseA,phaseB)
            comp = sys.getCompByName(compName);

            chems = struct2cell(sys.chemicals);
            if strcmp(phaseA,'Liquid')
                phaseACompSym = comp.getSym();
                phaseASym = sys.environs.(sys.activeEnv).getVComp().getSym();
            elseif strcmp(phaseA,'Gas')
                phaseACompSym = comp.getBulkGasSym();
                phaseASym = '(V_m-V_tot)';
            else
                for k=1:1:length(chems)
                    if strcmp(chems{k}.getType(),'Suspended Solid Sorbent') && strcmp(chems{k}.getName(),phaseA)
                        for l=1:1:length(comp.sorpFuncParams)
                            if strcmp(comp.sorpFuncParams{l}.solventName,phaseA)
                                phaseACompSym = comp.sorpFuncParams{l}.solventSym;
                            end
                        end
                        phaseASym = chems{k}.getSym();
                    end
                end
            end

            if strcmp(phaseB,'Liquid')
                phaseBSym = sys.environs.(sys.activeEnv).getVComp().getSym();
                if strcmp(phaseA,'Gas')
                    HConstHelper = comp.getHConstHelperFuncs("eq_with_liq");
                    phaseBCompSym = HConstHelper.getSubFuncSym();
                else
                    for k=1:1:length(chems)
                        if strcmp(chems{k}.getType(),'Suspended Solid Sorbent') && strcmp(chems{k}.getName(),phaseA)
                            for l=1:1:length(comp.sorpHelpers)
                                if strcmp(comp.sorpHelpers{l}.getSubFuncName(),[char(comp.name),' Concentration in ',phaseA,' Phase in Equilibrium with Liquid Phase'])
                                    phaseBCompSym = comp.sorpHelpers{l}.getSubFuncSym();
                                end
                            end
                        end
                    end
                end
            elseif strcmp(phaseB,'Gas')
                HConstHelper = comp.getHConstHelperFuncs("eq_with_gas");
                phaseBCompSym = HConstHelper.getSubFuncSym();
                phaseBSym = '(V_m-V_tot)';
            else
                for k=1:1:length(chems)
                    if strcmp(chems{k}.getType(),'Suspended Solid Sorbent') && strcmp(chems{k}.getName(),phaseB)
                        phaseBSym = chems{k}.getSym();
                        for l=1:1:length(comp.sorpHelpers)
                            if strcmp(comp.sorpHelpers{l}.getSubFuncName(),[char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',phaseB])
                                phaseBCompSym = comp.sorpHelpers{l}.getSubFuncSym();
                            end
                        end
                    end
                end
            end

            test = funcName
            test2 = phaseACompSym
            test3 = phaseBCompSym
            test4 = phaseASym
            test5 = phaseBSym
            MTFunc = CompDefaults.getDefaultMTFuncVals(funcName,phaseACompSym,phaseBCompSym,phaseASym,phaseBSym);
        end

        % sys: ODESys class ref, compName: string, phaseA: string, phaseB:
        % string
        function currModelType = getCurrMTFuncType(sys,compName,phaseA,phaseB)
            % ### STARTHERE: test this
            comp = sys.getCompByName(compName);
            if strcmp(phaseA,'Liquid')
                fp = comp.funcParams;
                for k=1:1:length(fp)
                    phaseBMatch = split(fp{k}.funcName," ");
                    if strcmp(phaseB,phaseBMatch{1})
                        currModelType = fp{k}.funcType;
                        return;
                    end
                end
            elseif strcmp(phaseA,'Gas')
                fp = comp.gasBulkFuncParams;
                for k=1:1:length(fp)
                    phaseBMatch = split(fp{k}.funcName," ");
                    if strcmp(phaseB,phaseBMatch{1})
                        currModelType = fp{k}.funcType;
                        return;
                    end
                end
            else
                fp = comp.sorpFuncParams;
                for k=1:1:length(fp)
                    phaseAMatch = split(fp{k}.funcName," ");
                    if strcmp(phaseA,phaseAMatch{1})
                        currModelType = fp{k}.funcType;
                        return;
                    end
                end
            end
            test1 = compName
            test2 = phaseA
            test3 = phaseB

        end

        % sys: ODESys class ref, strm_types: {}
        function sys = setStrmTypes(sys)
            sys.has_liq_in = false;
            sys.has_gas_in = false;
            sys.has_liq_out = false;
            sys.has_gas_out = false;
            for k=1:1:length(sys.input_streams)
                if strcmp(sys.input_streams{k}.getPhase(),'L'), sys.has_liq_in = true; end
                if strcmp(sys.input_streams{k}.getPhase(),'G'), sys.has_gas_in = true; end
            end
            for k=1:1:length(sys.output_streams)
                if strcmp(sys.output_streams{k}.getPhase(),'L'), sys.has_liq_out = true; end
                if strcmp(sys.output_streams{k}.getPhase(),'G'), sys.has_gas_out = true; end
            end
        end

        % sys: ODESys class ref, dir: string, phase: string
        function streams = getStrms(sys,dir,phase)
            streams = {};
            if dir == "in" && phase == "L"
                instreams = sys.input_streams;
                for k=1:1:length(instreams)
                    if instreams{k}.getPhase() == "L", streams{end+1} = instreams{k}; end %#ok<AGROW>
                end
            elseif dir == "in" && phase == "G"
                instreams = sys.input_streams;
                for k=1:1:length(instreams)
                    if instreams{k}.getPhase() == "G", streams{end+1} = instreams{k}; end %#ok<AGROW>
                end
            elseif dir == "out" && phase == "L"
                outstreams = sys.output_streams;
                for k=1:1:length(outstreams)
                    if outstreams{k}.getPhase() == "L", streams{end+1} = outstreams{k}; end %#ok<AGROW>
                end
            elseif dir == "out" && phase == "G"
                outstreams = sys.output_streams;
                for k=1:1:length(outstreams)
                    if outstreams{k}.getPhase() == "G", streams{end+1} = outstreams{k}; end %#ok<AGROW>
                end
            end
        end

        % sys: ODESys class ref, strmName: string
        function strm = getStrmByName(sys,strm)
            for k=1:1:length(sys.input_streams)
                if strcmp(strm,sys.input_streams{k}.getName())
                    strm = sys.input_streams{k};
                    return;
                end
            end
            for k=1:1:length(sys.output_streams)
                if strcmp(strm,sys.output_streams{k}.getName())
                    strm = sys.output_streams{k};
                    return;
                end
            end
        end

        % sys: ODESys classs ref
        function strmData = getStrmTData(sys)
            strmData = cell(size(sys.input_streams,1)+size(sys.output_streams,1),4);
            for k=1:1:length(sys.input_streams)
                strm = sys.input_streams{k};
                strmData{k,1} = strm.getName();
                strmData{k,2} = strm.getSym();
                strmData{k,3} = strm.getPhase();
                strmData{k,4} = strm.getDir();
            end
            for k=length(sys.input_streams)+1:1:length(sys.input_streams)+length(sys.output_streams)
                strm = sys.output_streams{k-length(sys.input_streams)};
                strmData{k,1} = strm.getName();
                strmData{k,2} = strm.getSym();
                strmData{k,3} = strm.getPhase();
                strmData{k,4} = strm.getDir();
            end
        end

        % sys: ODESys class ref, strmName: string, dir: string,
        % phase: string
        function addStrm(sys,strmName,dir,phase)
            strmNum = length(sys.input_streams) + length(sys.output_streams) + 1;
            if strcmp(dir,'in')
                sys.input_streams{end+1} = Stream(strmName,phase,dir,strmNum);
            elseif strcmp(dir,'out')
                sys.output_streams{end+1} = Stream(strmName,phase,dir,strmNum);
            end

            sys.inportCt = length(sys.input_streams);
            sys.outportCt = length(sys.output_streams);

            sys.setStrmTypes();
        end

        % sys: ODESys class ref, strmName: string, dir: string, phase:
        % string
        function removeStrm(sys,strmName,dir,phase)
            if dir == "in"
                for k=1:1:length(sys.input_streams)
                    if strcmp(sys.input_streams{k}.getName(),strmName) && ...
                            strcmp(sys.input_streams{k}.getPhase(),phase)
                        sys.input_streams(k) = [];
                        return;
                    end
                end
            elseif dir == "out"
                for k=1:1:length(sys.output_streams)
                    if strcmp(sys.output_streams{k}.getName(),strmName) && ...
                            strcmp(sys.output_streams{k}.getPhase(),phase)
                        sys.output_streams(k) = [];
                        return;
                    end
                end
            end

            sys.inportCt = length(sys.input_streams);
            sys.outportCt = length(sys.output_streams);

            sys.setStrmTypes();
        end

        % sys: ODESys class ref, name: string, phase: string, dir: string, compsData: {},
        % basis: string, basisU: string, solventName: string, flowrate: number, flowrateU: string
        function updateStrm(sys,name,dir,phase,compsData,basis,basisU,solventName,flowrate,flowrateU)
            if isempty(compsData)
                for k=1:1:length(sys.input_streams)
                    if strcmp(sys.input_streams{k}.getName(),name)
                        sys.input_streams{k}.setPhase(phase);
                        sys.input_streams{k}.setDir(dir);
                        return;
                    end
                end
                for k=1:1:length(sys.output_streams)
                    if strcmp(sys.output_streams{k}.getName(),name)
                        sys.output_streams{k}.setPhase(phase);
                        sys.output_streams{k}.setDir(dir);
                        return;
                    end
                end
            else
                for k=1:1:length(sys.input_streams)
                    if strcmp(sys.input_streams{k}.getName(),name)
                        sys.input_streams{k}.setCompsData(compsData,basis,basisU,solventName,flowrate,flowrateU);
                        return;
                    end
                end
                for k=1:1:length(sys.output_streams)
                    if strcmp(sys.output_streams{k}.getName(),name)
                        sys.output_streams{k}.setCompsData(compsData,basis,basisU,solventName,flowrate,flowrateU);
                        return;
                    end
                end
            end
        end

        % sys: ODESys class ref, strmName: string
        function strmCompsData = getStrmCompsData(sys,strmName)
            for k=1:1:length(sys.input_streams)
                if strcmp(sys.input_streams{k}.getName(),name)
                    strmCompsData = sys.input_streams{k}.getCompsData(strmName);
                    return;
                end
            end
            for k=1:1:length(sys.output_streams)
                if strcmp(sys.output_streams{k}.getName(),name)
                    strmCompsData = sys.output_streams{k}.getCompsData(strmName);
                    return;
                end
            end
        end

        % sys: ODESys class ref
        function names = getStrmSolvNames(sys)
            names = {'Water','V'};
            chems = struct2cell(sys.chemicals);
            if length(chems) > 0 %#ok<ISMT> 
                for k=1:1:length(chems)
                    if strcmp(chems{k}.getType(),'Suspended Solid Sorbent')
                        names{end+1,1} = char(chems{k}.getName()); %#ok<AGROW>
                        names{end,2} = char(chems{k}.getSym());
                    end
                end
            end
        end

        % sys: ODESys class ref
        % ### ACTUALLYSTARTHEREEXTENDED: need to differentiate between
        % liquid and suspended solid input streams
        function updateIOFunc(sys)
            % adds liquid volume IO to all relevant governing functions
            % V, X_j, C_j
            V_comp = sys.environs.(sys.activeEnv).getVComp();
            % ### FIXME: SV: for each SV, need to loop through and add the
            % proper terms for adjusting the volume due to I/O
            SV_comps = {};
            chems = sys.getChemicals('comp');
            for k=1:1:length(chems)
                if strcmp(chems{k}.getType(),'Suspended Solid Sorbent')
                    SV_comps{end+1} = chems{k}; %#ok<AGROW>
                end
            end
            % adding input/output functions for the SV
            % based on total liquid mass flow in and mass fraction of SV.
            % also need to create functions that modify the suspended solid
            % concentration, giving an average of the suspended solid
            % particles of different concentrations of solute
            % removes all I/O funcs from V_comp and SV_comps
            for k=2:1:length(V_comp.funcParams)
                fp = V_comp.funcParams{k};
                V_comp.removeModel(fp.funcVal,fp.funcName,'Main');
            end
            for k=1:1:length(SV_comps)
                for l=1:1:length(SV_comps{k}.funcParams)
                    fp = SV_comps{k}.funcParams{1};
                    SV_comps{l}.rmSuspendedSolidVolumeFunc(fp.funcVal,fp.funcName);
                end
            end
            % adds liquid I/O flow func variable syms to V_comp
            if ~sys.has_liq_in && ~sys.has_liq_out
                sys.addModel(V_comp.getName(),"0",'Liquid Flow Function','Main');
                for k=1:1:length(SV_comps)
                    sys.addSuspendedSolidVolumeFunc(SV_comps{k}.getName(),"0",'Liquid Flow Function');
                end
            else
                if sys.has_liq_in
                    sys.addModel(V_comp.getName(),"L_i",'Liquid Flow In','Main');
                    for k=1:1:length(SV_comps)
                        sys.addSuspendedSolidVolumeFunc(SV_comps{k}.getName(),"L_s"+SV_comps{k}.getNum()+"_i",'Liquid Flow In');
                    end
                end
                if sys.has_liq_out
                    sys.addModel(V_comp.getName(),"-L_o",'Liquid Flow Out','Main');
                    for k=1:1:length(SV_comps)
                        sys.addSuspendedSolidVolumeFunc(SV_comps{k}.getName(),"-L_s"+SV_comps{k}.getNum()+"_o",'Liquid Flow Out');
                    end
                end
            end
            % loops through each liquid-phase component, checks whether
            % liquid flow I/O func variable syms need to be added or
            % removed
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'), ...
                {sys.environs.(sys.activeEnv).getH3OComp(),sys.environs.(sys.activeEnv).getOHComp()}];
            for k=1:1:length(comps)
                fp2 = comps{k}.funcParams;
                has_lfi = false;
                for l=1:1:length(fp2)
                    if strcmp(fp2{l}.funcName,'Liquid Flow In')
                       has_lfi = true;
                    end
                end
                % adding term for liquid input flow into each component
                if sys.has_liq_in
                    if ~has_lfi
                        sys.addModel(comps{k}.getName(),"L_i/V*"+comps{k}.getSym()+"_i",'Liquid Flow In','Liquid');
                        for l=1:1:length(SV_comps)
                            if any(strcmp(SV_comps{l}.sorpTInfo(:,1),comps{k}.getName()))
                                sys.addModel(comps{k}.getName(),"L_s"+SV_comps{l}.getNum()+"_i/V_s*"+comps{k}.getSym()+"_i",'Liquid Flow In','Suspended Solid');
                            end
                        end
                    end
                else
                    % removing term for liquid input flow from each
                    % component
                    if has_lfi
                        comps{k}.removeModel("L_i/V*"+comps{k}.getSym()+"_i",'Liquid Flow In','Liquid');
                        for l=1:1:length(SV_comps)
                            if any(strcmp(SV_comps{l}.sorpTInfo(:,1),comps{k}.getName()))
                                comps{k}.removeModel("L_s"+SV_comps{l}.getNum()+"_i/V_s*S_"+comps{k}.getNum()+ ...
                                    "_"+SV_comps{l}.getNum()+"_i",'Liquid Flow In','Suspended Solid');
                            end
                        end
                    end
                end
                has_lfo = false;
                for l=1:1:length(fp2)
                    if strcmp(fp2{l}.funcName,'Liquid Flow Out')
                       has_lfo = true;
                    end
                end
                % adding term for liquid output flow to each component
                if sys.has_liq_out
                    if ~has_lfo
                        sys.addModel(comps{k}.getName(),"-L_o/V*"+comps{k}.getSym(),'Liquid Flow Out','Liquid');
                        for l=1:1:length(SV_comps)
                            if any(strcmp(SV_comps{l}.sorpTInfo(:,1),comps{k}.getName()))
                                sys.addModel(comps{k}.getName(),"-L_s"+SV_comps{l}.getNum()+"_o/V_s*S_"+comps{k}.getNum()+ ...
                                    "_"+SV_comps{l}.getNum(),'Liquid Flow In','Suspended Solid');
                            end
                        end
                    end
                else
                    % removing term for liquid output flow from each
                    % component
                    if has_lfo
                        comps{k}.removeModel("-L_o/V*"+comps{k}.getSym(),'Liquid Flow Out','Liquid');
                        for l=1:1:length(SV_comps)
                            if any(strcmp(SV_comps{l}.sorpTInfo(:,1),comps{k}.getName()))
                                comps{k}.removeModel("-L_s"+SV_comps{l}.getNum()+"_o/V_s*S_"+comps{k}.getNum()+ ...
                                    "_"+SV_comps{l}.getNum(),'Liquid Flow In','Suspended Solid');
                            end
                        end
                    end
                end
            end

            % adds gas volume IO to all relevant governing functions
            % P, P_j
            P_comp = sys.environs.(sys.activeEnv).getPComp();
            chems = sys.getChemicals('comp');
            for k=1:1:length(chems)
                fp3 = chems{k}.funcParams;
                has_gfi = false;
                for l=1:1:length(fp3)
                    if strcmp(fp3{l}.funcName,'Gas Flow In')
                       has_gfi = true;
                    end
                end
                if sys.has_gas_in
                    if ~has_gfi
                        % ### UPDATE: need to be updating the model in the
                        % chemicals and P_comp
                        if chems{k}.is_vol
                            % adding function to P_comp for bub flow input
                            % to bulk gas phase
                            sys.addModel(P_comp.getName(), ...
                            "(R*T)/((V_m-V_tot)*"+chems{k}.MW_sym+")*(D_i_+"+chems{k}.getSym()+"-(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Main');
                            % adding function to chem for bub flow input
                            % to bulk gas phase
                            sys.addModel(chems{k}.getName(), ...
                            "(R*T)/((V_m-V_tot)*"+chems{k}.MW_sym+")*(D_i_+"+chems{k}.getSym()+"-(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Gas');
                            % adding helper function for final bub
                            % concentration
                            C_bub_f_val = chems{k}.sym+"+(P_i*Y_i_"+chems{k}.number+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+")*exp(-tau_"+chems{k}.sym+"*k_"+chems{k}.sym+")";
                            sys.addRmHelperFuncs([char(chems{k}.name),' Bubble Final Concentration'],char(chems{k}.sym+"_bubf"),C_bub_f_val,true);
                            % adding function to chem for bub to liq
                            sys.addModel(chems{k}.getName(), ...
                            "(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Liquid');
                        end
                    end
                else
                    if has_gfi
                        if chems{k}.is_vol
                            % removing function from P_comp for bub flow input
                            % to bulk gas phase
                            P_comp.removeModel("(R*T)/((V_m-V_tot)*"+chems{k}.MW_sym+")*(D_i_+"+chems{k}.getSym()+"-(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Main');
                            % removing function from chem for bub flow input
                            % to bulk gas phase
                            chems{k}.removeModel("(R*T)/((V_m-V_tot)*"+chems{k}.MW_sym+")*(D_i_+"+chems{k}.getSym()+"-(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Gas');
                            % removing helper function for final bub
                            % concentration
                            C_bub_f_val = chems{k}.sym+"+(P_i*Y_i_"+chems{k}.number+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+")*exp(-tau_"+chems{k}.sym+"*k_"+chems{k}.sym+")";
                            sys.addRmHelperFuncs([char(chems{k}.name),' Bubble Final Concentration'],char(chems{k}.sym+"_bubf"),C_bub_f_val,false);
                            % removing function from chem for bub to liq
                            chems{k}.removeModel("(P_i*Y_i_"+chems{k}.sym+"/"+chems{k}.h_const_sym+"-"+chems{k}.sym+"_bubf))/tau*V", ...
                                'Gas Flow In','Liquid');
                        end
                    end
                end
                has_gfo = false;
                for l=1:1:length(fp3)
                    if strcmp(fp3{l}.funcName,'Gas Flow Out')
                       has_gfo = true;
                    end
                end
                if sys.has_gas_out
                    if ~has_gfo
                        if chems{k}.is_vol
                            % adding function to P_comp for bub flow output
                            % from bulk gas phase
                            sys.updateModel(P_comp.getName(),"-D_o*R*T/(V_m-V_tot)*"+chems{k}.bulk_gas_sym+"/(P*"+chems{k}.MW_sym+")", ...
                                'Gas Flow Out','Main');
                            % adding function to chem for bub flow output
                            % from bulk gas phase
                            sys.updateModel(chems{k}.getName(),"-D_o*R*T/(V_m-V_tot)*"+chems{k}.bulk_gas_sym+"/(P*"+chems{k}.MW_sym+")", ...
                                'Gas Flow Out','Gas');
                        end
                    end
                else
                    if has_gfo
                        % removing function from P_comp for bub flow output
                        % to bulk gas phase
                        P_comp.removeModel("-D_o*R*T/(V_m-V_tot)*"+chems{k}.bulk_gas_sym+"/(P*"+chems{k}.MW_sym+")",'Gas Flow Out','Main');
                        % removing function from chem for bub flow output
                        % to bulk gas phase 
                        chems{k}.removeModel("-D_o*R*T/(V_m-V_tot)*"+chems{k}.bulk_gas_sym+"/(P*"+chems{k}.MW_sym+")",'Gas Flow Out','Gas');
                    end
                end
            end
        end

        % sys: ODESys class ref
        function updateLiqVolEvts(sys)
            % ### STARTHERE: will need to run this whenever input/output
            % streams are added/removed, and whenever the flowrates for
            % streams are modified
        end

        % sys: ODESys class ref, solventName: string, soluteName: string,
        % modelName: string, param1Val: num, param2Val: num, param3Val: num
        function sorpTInfo = addSorptionEq(sys,solventName,soluteName,modelName,param1Val,param2Val,param3Val)
            solvComp = sys.getCompByName(solventName);
            soluComp = sys.getCompByName(soluteName);
            solventNum = sys.getChemicalNumFromName(solventName);
            soluteNum = sys.getChemicalNumFromName(soluteName);
            % keeps track of governing functions for equilibrium and mass
            % transfer
            sorpHelpers = soluComp.addSorptionEq(solventName,solventNum,soluteName,soluteNum,modelName,param1Val,param2Val,param3Val,sys.getDefaultParamVals());
            % keeps track of sorpTInfo
            sorpTInfo = solvComp.addSorptionEq(soluteName,solventNum,soluteName,soluteNum,modelName,param1Val,param2Val,param3Val,sys.getDefaultParamVals());

            % add sorption helpers
            helperFuncNames = sys.getHelperFuncNames();
            for k=1:1:length(sorpHelpers)
                if ~any(strcmp(helperFuncNames,sorpHelpers{k}.getSubFuncName()))
                    sys.addRmHelperFuncs(sorpHelpers{k}.funcName,sorpHelpers{k}.funcSym,sorpHelpers{k}.funcVal,true);
                end
                [~,~,helperIdx] = sys.getHelperFunc(sorpHelpers{k}.getSubFuncName());
                paramNames = sorpHelpers{k}.getSubFuncParamNames();
                paramSyms = sorpHelpers{k}.getSubFuncParamSyms();
                paramVals = sorpHelpers{k}.getSubFuncParamVals();
                paramUnits = sorpHelpers{k}.getSubFuncParamUnits();
                for l=1:1:length(paramNames)
                    sys.helperFuncs{helperIdx}.updateParams(paramNames(l),paramSyms(l),paramVals(l),paramUnits(l),true,sys.getDefaultParamVals());
                end
            end
            % test4 = sys.helperFuncs
        end

        % sys: ODESys class ref, solventName: string, soluteName: string,
        % modelName: string, param1Val: num, param2Val: num, param3Val: num
        function sorpTInfo = rmSorptionEq(sys,solventName,soluteName,modelName,param1Val,param2Val,param3Val)
            solvComp = sys.getCompByName(solventName);
            soluComp = sys.getCompByName(soluteName);
            solventNum = sys.getChemicalNumFromName(solventName);
            % returns helper functions that need to be removed
            sorpHelpers = soluComp.rmSorptionEq(solventName,solventNum,soluteName,modelName,param1Val,param2Val,param3Val);
            % soluteNum = sys.getChemicalNumFromName(soluteName);
            % keeps track of sorpTInfo
            sorpTInfo = solvComp.rmSorptionEq(solventName,solventNum,soluteName,modelName,param1Val,param2Val,param3Val);

            % remove sorption helpers
            for k=1:1:length(sorpHelpers)
                sys.addRmHelperFuncs(sorpHelpers{k}.funcName,sorpHelpers{k}.funcSym,sorpHelpers{k}.funcVal,false);
            end
        end

        % sys: ODESys class ref
        function res = getInteractCompDDItems(sys)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),{sys.environs.(sys.activeEnv).getTComp(), ...
                sys.environs.(sys.activeEnv).getH3OComp(),sys.environs.(sys.activeEnv).getOHComp()}];
            res = cell(length(comps),2);
            for k=1:1:length(comps)
                if length(comps) - k >= 5
                    name = [char(comps{k}.getName()),' (Liquid Phase)'];
                else
                    name = char(comps{k}.getName());
                end
                res{k,1} = name;
                res{k,2} = comps{k}.getSym();
            end
        end

        % sys: ODESys class ref
        function res = getCompParamDDItems(sys)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp')];
            res = cell(length(comps)+4,2);
            res{1,1} = 'Temperature';
            res{1,2} = 'T';
            res{2,1} = 'Pressure';
            res{2,2} = 'P';
            res{3,1} = 'Volume';
            res{3,2} = 'V';
            res{4,1} = 'Hydronium';
            res{4,2} = 'H3O';
            res{5,1} = 'Hydroxide';
            res{5,2} = 'OH';
            for k=1:1:length(comps)
                res{k+5,1} = char(comps{k}.getName());
                res{k+5,2} = comps{k}.getSym();
            end
            for k=1:1:length(sys.helperFuncs)
                res{k+5+length(comps),1} = char(sys.helperFuncs{k}.getSubFuncName());
                res{k+5+length(comps),2} = char(sys.helperFuncs{k}.getSubFuncSym());
            end
        end

        % sys: ODESys class ref
        function res = getRegFuncNames(sys)
            comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
            res = cell(length(comps),2);
            for k=1:1:length(comps)
                if any(contains({'T','P','V'},comps{k}.getSym()))
                    ext = '';
                else
                    ext = ' (Liquid Phase)';
                end
                res{k,1} = [char(comps{k}.getName()),ext];
                res{k,2} = comps{k}.getSym();
            end

            subfuncs = [sys.environs.(sys.activeEnv).subfuncs,sys.helperFuncs];
            for k=1:1:length(subfuncs)
                res{k+length(comps),1} = char(subfuncs{k}.getSubFuncName());
                res{k+length(comps),2} = subfuncs{k}.getSubFuncSym();
            end
        end

        % sys: ODESys class ref
        function runModelAndExportData(sys)
            % compile model
            sys.compileModel();

            % generate data based on existing plots
            tEnd = sys.model_runtime;
            tSmooth = linspace(0,tEnd,100);
            sysVar = sys.getModelVarNames("plot");
            for k=1:1:length(sys.plots)
                plot_obj = sys.plots{k};
                axes = plot_obj.axes;
                if length(axes) == 3
                    if axes{1}.varIsIC && axes{2}.varIsIC
                        y0_span_x = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts);
                        y0_span_y = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts);
                        % maybe remove these and just use y0_mesh_X and
                        % y0_mesh_Y for the interpolation query points
                        % as well?
                        y0_span_qx = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts*5);
                        y0_span_qy = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts*5);
                        [y0_mesh_X,y0_mesh_Y] = meshgrid(y0_span_x,y0_span_y);
                        % [y0_mesh_qX,y0_mesh_qY] = meshgrid(y0_span_qx,y0_span_qy);

                        res_1 = cell(length(y0_span_y),length(y0_span_x));

                        % creating array that maps y0 to variable name
                        y0_2_name = [];
                        y0_other_phase = [];
                        comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                        ext_ct = 0;
                        for l=1:1:length(comps)
                            y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Liquid Phase)']); %#ok<AGROW>
                            if comps{l}.is_vol
                                ext_ct = ext_ct + 1;
                                % helperName = ['Gas Phase Partial Pressure of ',char(comps{l}.getName()),' in Equilibrium with Liquid Phase'];
                                % [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                % eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                % % ### FIXME?: unable to access initial conditions for
                                % % variables inputted from Simulink, will assume 0 for
                                % % all
                                % sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));

                                y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Gas Phase)']); %#ok<AGROW>
                                y0_other_phase(end+1) = l; %#ok<AGROW>
                            end
                            solvents = {};
                            if ~isempty(comps{k}.sorpFuncParams)
                                for m=1:1:length(comps{k}.sorpFuncParams)
                                    if ~any(strcmp(solvents,comps{k}.sorpFuncParams{l}.solventName))
                                        ext_ct = ext_ct + 1;
                                        % helperName = [char(comps{k}.name),' Concentration in ',comps{k}.sorpFuncParams{l}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                        % [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        % eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        % sys.y0(k+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));

                                        y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (',comps{k}.sorpFuncParams{l}.solventName,' Phase)']); %#ok<AGROW>
                                        y0_other_phase(end+1) = l; %#ok<AGROW>
                                        solvents{end+1} = comps{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                                    end
                                end
                            end
                        end

                        y0_var_num_x = strcmp(y0_2_name,erase(axes{1}.getPlotProp("varNames"),"~(IC)"));
                        y0_var_num_y = strcmp(y0_2_name,erase(axes{2}.getPlotProp("varNames"),"~(IC)"));

                        multiple_x = length(find(y0_other_phase(y0_var_num_x) == y0_other_phase));
                        multiple_y = length(find(y0_other_phase(y0_var_num_y) == y0_other_phase));

                        for l=1:1:length(y0_span_y)
                            for m=1:1:length(y0_span_x)
                                y0_1 = sys.y0;
                                % ### FIXME: need to write logic to edit the
                                % changes in initial conditions when there are
                                % chemical solutes in different phases
                                y0_1(y0_var_num_x) = y0_span_(m);
                                if multiple_x
                                    first_x_idx = find(y0_other_phase(y0_var_num_x) == y0_other_phase,1);
                                    % setting comp
                                    comp = comps{y0_other_phase(first_x_idx)};
                                    if first_x_idx + 1 == y0_var_num_x
                                        helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_x_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    elseif first_x_idx ~= y0_var_num_x
                                        for n=1:1:length(comp.sorpFuncParams)
                                            if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{n}.solventName)
                                                solventName = comp.sorpFuncParams{n}.solventName;
                                            end
                                        end
                                        helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_x_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end

                                    % now update the other phases
                                    ext_ct = 0;
                                    if comp.is_vol
                                        ext_ct = ext_ct + 1;
                                        helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        % ### FIXME?: unable to access initial conditions for
                                        % variables inputted from Simulink, will assume 0 for
                                        % all
                                        y0_1(first_x_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end
                                    solvents = {};
                                    if ~isempty(comp.sorpFuncParams)
                                        for n=1:1:length(comp.sorpFuncParams)
                                            if ~any(strcmp(solvents,comp.sorpFuncParams{n}.solventName))
                                                ext_ct = ext_ct + 1;
                                                helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{n}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                                [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                                eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                                y0_1(first_x_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                                solvents{end+1} = comp.sorpFuncParams{n}.solventName; %#ok<AGROW>
                                            end
                                        end
                                    end
                                end
                                y0_1(y0_var_num_y) = y0_span_y(l);
                                if multiple_y
                                    first_y_idx = find(y0_other_phase(y0_var_num_y) == y0_other_phase,1);
                                    % setting comp
                                    comp = comps{y0_other_phase(first_y_idx)};
                                    if first_y_idx + 1 == y0_var_num_y
                                        helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_y_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    elseif first_y_idx ~= y0_var_num_y
                                        for n=1:1:length(comp.sorpFuncParams)
                                            if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{n}.solventName)
                                                solventName = comp.sorpFuncParams{n}.solventName;
                                            end
                                        end
                                        helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        y0_1(first_y_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end

                                    % now update the other phases
                                    ext_ct = 0;
                                    if comp.is_vol
                                        ext_ct = ext_ct + 1;
                                        helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                        [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                        eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                        % ### FIXME?: unable to access initial conditions for
                                        % variables inputted from Simulink, will assume 0 for
                                        % all
                                        y0_1(first_y_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                    end
                                    solvents = {};
                                    if ~isempty(comp.sorpFuncParams)
                                        for n=1:1:length(comp.sorpFuncParams)
                                            if ~any(strcmp(solvents,comp.sorpFuncParams{n}.solventName))
                                                ext_ct = ext_ct + 1;
                                                helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{n}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                                [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                                eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                                y0_1(first_y_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                                solvents{end+1} = comp.sorpFuncParams{n}.solventName; %#ok<AGROW>
                                            end
                                        end
                                    end
                                end

                                [tRes,yRes] = sys.runModel(tSmooth,y0_1);
                                fRes = sys.calculateHelperVals(tRes,yRes);
                                res_1{l,m} = [yRes,fRes,tRes];
                            end
                        end

                        res_2 = cell(length(axes{3}.getPlotProp("varNames")),1);
                        evalt = sys.plots{k}.getPlotProp("evaltVal");
                        for l=1:1:length(axes{3}.getPlotProp("varNames"))
                            yVarIdx = strcmp(sysVar(:,1),axes{3}.varNames{l});
                            res_2{l} = zeros(size(res_1));
                            for m=1:1:size(res_2,1)
                                for n=1:1:size(res_2,2)
                                    % ### FIXME: add more methods for
                                    % interpolation?
                                    res_2{l}(m,n) = interp1(res_1{m,n}(:,end),res_1{m,n}(:,yVarIdx),evalt,'makima');
                                end
                            end
                        end

                        % ### FIXME: need to properly do this
                        dataT = table(y0_span_x,y0_span_y,'VariableNames',{axes{1}.getPlotProp("title"),axes{2}.getPlotProp("title")});
                        for l=1:1:length(res_2)
                            dataT = [dataT,table(y0_span_x,res_2{l},'VariableNames',{['Columns: ',axes{2}.getPlotProp("title"),'/Rows: ',axes{1}.getPlotProp("title")], ...
                                y0_span_y})]; %#ok<AGROW>
                        end
                    elseif axes{1}.varIsIC || axes{2}.varIsIC
                        if axes{1}.varIsIC
                            y0_span = linspace(axes{1}.loEvalLim,axes{1}.upEvalLim,axes{1}.nbEvalPts);
                            ICax = 1;
                            ICax_obj = axes{1};
                            nonICAx = 2;
                            nonICAx_obj = axes{2};
                        else
                            y0_span = linspace(axes{2}.loEvalLim,axes{2}.upEvalLim,axes{2}.nbEvalPts);
                            ICax = 2;
                            ICax_obj = axes{2};
                            nonICAx = 1;
                            nonICAx_obj = axes{1};
                        end

                        res_1 = cell(length(y0_span),1);

                        % creating array that maps y0 to variable name
                        y0_2_name = [];
                        y0_other_phase = [];
                        comps = [sys.getSpecies('comp'),sys.getChemicals('comp'),sys.environs.(sys.activeEnv).getAllEnvComps()];
                        ext_ct = 0;
                        for l=1:1:length(comps)
                            y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Liquid Phase)']); %#ok<AGROW>
                            if comps{l}.is_vol
                                ext_ct = ext_ct + 1;
                                y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (Gas Phase)']); %#ok<AGROW>
                                y0_other_phase(end+1) = l; %#ok<AGROW>
                            end
                            solvents = {};
                            if ~isempty(comps{k}.sorpFuncParams)
                                for m=1:1:length(comps{k}.sorpFuncParams)
                                    if ~any(strcmp(solvents,comps{k}.sorpFuncParams{l}.solventName))
                                        ext_ct = ext_ct + 1;
                                        y0_2_name(l+ext_ct) = string([comps{l}.getName(),' (',comps{k}.sorpFuncParams{l}.solventName,' Phase)']); %#ok<AGROW>
                                        y0_other_phase(end+1) = l; %#ok<AGROW>
                                        solvents{end+1} = comps{k}.sorpFuncParams{l}.solventName; %#ok<AGROW>
                                    end
                                end
                            end
                        end

                        y0_var_num = strcmp(y0_2_name,erase(plot_obj.getAxProp(ICax,"varNames"),"~(IC)"));

                        multiple_num = length(find(y0_other_phase(y0_var_num) == y0_other_phase));

                        for l=1:1:length(y0_span)
                            y0_1 = sys.y0;
                            y0_1(y0_var_num) = y0_span(l);
                            if multiple_num
                                first_idx = find(y0_other_phase(y0_var_num) == y0_other_phase,1);
                                % setting comp
                                comp = comps{y0_other_phase(first_idx)};
                                if first_idx + 1 == y0_var_num
                                    helperName = ['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'];
                                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                    y0_1(first_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                elseif first_idx ~= y0_var_num
                                    for m=1:1:length(comp.sorpFuncParams)
                                        if contains(axes{1}.getPlotProp("varNames"),comp.sorpFuncParams{m}.solventName)
                                            solventName = comp.sorpFuncParams{m}.solventName;
                                        end
                                    end
                                    helperName = [char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName];
                                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                    y0_1(first_idx) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                end

                                % now update the other phases
                                ext_ct = 0;
                                if comp.is_vol
                                    ext_ct = ext_ct + 1;
                                    helperName = ['Gas Phase Partial Pressure of ',char(comp.getName()),' in Equilibrium with Liquid Phase'];
                                    [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                    eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                    y0_1(first_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                end
                                solvents = {};
                                if ~isempty(comp.sorpFuncParams)
                                    for m=1:1:length(comp.sorpFuncParams)
                                        if ~any(strcmp(solvents,comp.sorpFuncParams{m}.solventName))
                                            ext_ct = ext_ct + 1;
                                            helperName = [char(comp.name),' Concentration in ',comp.sorpFuncParams{m}.solventName,' Phase in Equilibrium with Liquid Phase'];
                                            [~,~,helperIdx] = sys.getHelperFunc(helperName);
                                            eqFunc = sys.f{helperIdx+length(sys.environs.(sys.activeEnv).subfuncs),3};
                                            y0_1(first_idx+ext_ct) = eqFunc(0,sys.y0,sys.param,sys.f(:,3),zeros(size(sys.var_in_key)));
                                            solvents{end+1} = comp.sorpFuncParams{m}.solventName; %#ok<AGROW>
                                        end
                                    end
                                end

                                [tRes,yRes] = sys.runModel(tSmooth,y0_1);
                                fRes = sys.calculateHelperVals(tRes,yRes);
                                res_1{l,1} = [yRes,fRes,tRes];
                            end
                        end

                        res_2 = {};
                        nonICVarIdx = strcmp(sysVar(:,1),nonICAx_obj.varNames);
                        for l=1:1:length(axes{3}.getPlotProp("varNames"))
                            zVarIdx = strcmp(sysVar(:,1),axes{3}.varNames{l});
                            res_2{l} = zeros(size(res_1)); %#ok<AGROW>
                            for m=1:1:size(res_2,1)
                                % ### FIXME: add more methods for
                                % interpolation?
                                res_2{l}(m,1) = interp1(res_1{m,1}(:,nonICVarIdx),res_1{m,1}(:,zVarIdx), ...
                                    linspace(min(res_1{m,1}(:,nonICVarIdx)),max((res_1{m,1}(:,nonICVarIdx)),size(res_1{m,1},1))),'makima');
                            end
                        end

                        % ### DEBUG: need to try this
                        dataT = table(y0_span_x,y0_span_y,'VariableNames',{axes{1}.getPlotProp("title"),axes{2}.getPlotProp("title")});
                        for l=1:1:length(res_2)
                            dataT = [dataT,table(y0_span_x,res_2{l},'VariableNames',{['Columns: ',axes{2}.getPlotProp("title"),'/Rows: ',axes{1}.getPlotProp("title")], ...
                                y0_span_y})]; %#ok<AGROW>
                        end
                    else
                        [tRes,yRes] = sys.runModel(tSmooth);
                        fRes = sys.calculateHelperVals(tRes,yRes);
                        res_1{1,1} = [yRes,fRes,tRes];

                        xVarIdx = strcmp(sysVar(:,1),axes{1}.varNames);
                        yVarIdx = strcmp(sysVar(:,1),axes{2}.varNames);
                        zVarIdx = zeros(size(axes{3}.varNames));
                        for l=1:1:length(axes{3}.varNames)
                            zVarIdx(l) = strcmp(sysVar(:,1),axes{2}.varNames{l});
                        end

                        y0_span_x = res_1{1,1}(:,xVarIdx);
                        y0_span_y = res_1{1,1}(:,yVarIdx);
                        y0_span_z = res_1{1,1}(:,zVarIdx);

                        % creating data export table
                        dataT = table(y0_span_x,y0_span_y,y0_span_z,'VariableNames',{sysVar(xVarIdx,1),sysVar(yVarIdx,1),sysVar(zVarIdx,1)});
                    end
                else
                    [tRes,yRes] = sys.runModel(tSmooth);
                    fRes = sys.calculateHelperVals(tRes,yRes);
                    res_1{1,1} = [yRes,fRes,tRes];

                    % ### FIXME: export data with different time units
                    dataT = table(res_1{1,1},'VariableNames',{sysVar(end,1),sysVar(1:end-1,1)});
                end

                % write table
                writetable(dataT,sys.data_export_dir,'Sheet',sys.plots{k}.getPlotProp('title'));
            end
        end

        % sys: ODESys class ref
        function solver = getODESolver(sys)
            solver = sys.solver;
        end

        % sys: ODESys class ref, newSolver: ode
        function setODESolver(sys,newSolver)
            sys.solver = newSolver;
        end
    end

    methods (Static)
        % funcAsStr: Boolean
        function [paramVals, paramUnits] = getDefaultEnvironParams(funcAsStr)
            custDefault = EnvDefaults();
            paramVals = custDefault.Values;
            paramUnits = custDefault.Units;
            if funcAsStr
                for i=1:1:length(paramVals)
                    if isa(paramVals{i},"function_handle")
                        paramVals{i} = func2str(paramVals{i});
                    end
                end
            end
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