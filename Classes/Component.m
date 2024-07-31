classdef Component < handle
    properties
        number = 0; % comps identifying number
        name = ''; % comps name

        initConc = 1; % initial concentration for comps
        initConcUnit = '';

        funcParams = {}; % cell array of structs that contain the values, names, symbols, and units of parameters,
                        % with the indicies being the parameter numbers.
        sym = '';
        type = "";

        % Liquid-gas mass transfer properties
        h_const = 1; % Henry's Law Constant
        h_const_u = ''; % Henry's Law Constant unit
        h_const_sym = ''; % Henry's Law Constant symbol
        dh_const = 1; % Temperature Dependence Coefficient for Henry's Law
        dh_const_u = '';
        dh_const_sym = '';
        h_const_helper_funcs = {};
        is_vol = false;
        bulk_gas_sym = "";
        gasBulkFuncParams = {};
        gasBulkHelpers = {};

        % Suspended-solid phase mass transfer properties
        sorped_sym = struct([]);
        % sorpFuncParams is same structure as funcParams, but has extra
        % properties: solventName: string, solventSym: string
        % in compileGovFunc, need to loop through every sorpFuncParam to
        % create a new dependent variable for each phase
        sorpFuncParams = {};
        sorpVolFuncParams = {};
        sorpHelpers = {};
        sorpTInfo = {};

        MW = 1; % Molecular Weight
        MW_sym = '';
        den = 1; % Density
        den_u = ''; % Density unit
        den_sym = '';
        Cp = 1; % Heat Capacity
        Cp_u = ''; % Heat Capacity unit
        Cp_sym = '';
        Cpg = 1; % Gas Phase Heat Capacity
        Cpg_u = ''; % Gas Phase Heat Capacity unit
        Cpg_sym = '';

        %% funcParams Item Data Structure: {}
        %   funcVal: string
        %   funcName: string
        %   params: {struct}
        %       locNum: number
        %       gloNum: number
        %       name: string
        %       sym: string
        %       val: number
        %       unit: string
        %       editable: boolean
        %   lims: {struct - by variable name}
        %       funcVal: string
        %       lowerLim: number
        %       upperLim: number
        %       sysVar: string
        %       elseVal: number
        %       lowLimOp: string
        %       upLimOp: string

    end

    methods
        % varargin:
        % 1 - MW
        % 2 - den
        % 3 - den_u
        % 4 - Cp
        % 5 - Cp_u
        % 6 - h_const
        % 7 - h_const_u
        % 8 - dh_const
        % 9 - dh_const_u
        function comp = Component(name,number,initConc,initConcUnit,type,is_vol,MW,h_const,h_const_u,dh_const,dh_const_u,defaultParamVals)
            comp.name = name;
            comp.number = number;
            if strcmp(type,'Biological Solute')
                comp.sym = char("X_"+number);
                % creates funcParam for biological governing functions
                comp.funcParams{1} = comp.createFuncParamObj("0",'BioGovFunc');
            elseif strcmp(type,'Chemical Solute')
                comp.sym = char("C_"+number);
                % creates funcParam for biological governing functions
                comp.funcParams{1} = comp.createFuncParamObj("0",'BioGovFunc');
            elseif strcmp(type,'Suspended Solid Sorbent')
                comp.sym = char("V_S_"+number);
                comp.funcParams{1} = comp.createFuncParamObj("0",'MTFunc');
            elseif strcmp(type,'Liquid Solvent')
                comp.sym = char("LS_"+number);
                comp.funcParams{1} = comp.createFuncParamObj("0",'MTFunc');
            elseif strcmp(type,'temp')
                comp.sym = 'T';
                comp.funcParams{1} = comp.createFuncParamObj("0",'TempFunc');
            elseif strcmp(type,'press')
                comp.sym = 'P';
                comp.funcParams{1} = comp.createFuncParamObj("0",'PressFunc');
            elseif strcmp(type,'vol')
                comp.sym = 'V';
                comp.funcParams{1} = comp.createFuncParamObj("0",'VolFunc');
            elseif strcmp(type,'acid')
                comp.sym = 'H3O';
                initConc = initConc * MW;
                comp.funcParams{1} = comp.createFuncParamObj("0",'AcidFunc');
            elseif strcmp(type,'base')
                comp.sym = 'OH';
                initConc = initConc * MW;
                comp.funcParams{1} = comp.createFuncParamObj("0",'BaseFunc');
            end
            comp.initConc = initConc;
            comp.initConcUnit = initConcUnit;
            comp.type = type;

            % set species parameters
            comp.MW = MW;
            comp.MW_sym = ['MW_',char(comp.sym)];

            % comp.den_u = den_u;
            % switch comp.den_u
            %     case 'g/L'
            %         comp.den = den;
            %     case 'kg/m^3'
            %         comp.den = den;
            %     case 'lbm/ft^3'
            %         comp.den = den / 28.168 * 453.592;
            % end
            % comp.den_sym = ['rho_',char(comp.sym)];
            % 
            % comp.Cp_u = Cp_u;
            % comp.Cpg_u = Cp2_u;
            % switch comp.Cp_u
            %     case 'kJ/(kg*K)'
            %         comp.Cp = Cp;
            %         comp.Cpg = Cp2;
            %     case 'J/(g*K)'
            %         comp.Cp = Cp;
            %         comp.Cpg = Cp2;
            % end
            % comp.Cp_sym = ['C_p_',char(comp.sym)];
            % comp.Cpg_sym = ['C_pg_',char(comp.sym)];

            if is_vol, comp.setVol(is_vol,h_const,h_const_u,dh_const,dh_const_u,defaultParamVals); end
            comp.is_vol = is_vol;
        end

        % comp: Compnent class ref, depVars: char{}
        function comp = setHConstHelperFuncs(comp,defaultParamVals)
            % sets helpers for the Henry's Constant function, vapor
            % partial pressure in equilibrium with liquid phase
            % concentration, and liquid phase concentration in equilibrium
            % with vapor partial pressure
            % Total 3 helpers (HConst, P*, C*)
            % Need to consider water density?
            funcVal = [comp.h_const_sym,'*',comp.MW_sym,'*exp(',comp.dh_const_sym,'*((1/T)-(1/298.15)))'];
            comp.h_const_helper_funcs{end+1} = SubFunc(funcVal,[char(comp.name),' Henry Constant'],['H_',char(comp.sym)],0,0);
            comp.h_const_helper_funcs{end}.updateParams([char(comp.name),' Henry Constant at 298.15K'],comp.h_const_sym,comp.h_const,comp.h_const_u,false,defaultParamVals);
            comp.h_const_helper_funcs{end}.updateParams([char(comp.name),' H Temperature Dependence Coefficient'],comp.dh_const_sym,comp.dh_const,comp.dh_const_u,false,defaultParamVals);

            funcVal2 = [char(comp.sym),'/H_',char(comp.sym)];
            comp.h_const_helper_funcs{end+1} = SubFunc(funcVal2,['Gas Phase Partial Pressure of ',char(comp.name),' in Equilibrium with Liquid Phase'], ...
                ['P_eq_',char(comp.sym)],0,0);
            
            funcVal3 = [comp.bulk_gas_sym,'*H_',char(comp.sym)];
            comp.h_const_helper_funcs{end+1} = SubFunc(funcVal3,['Liquid Phase Concentration of ',char(comp.name),' in Equilibrium with Gas Phase'], ...
                [char(comp.sym),'_g_eq'],0,0);
        end

        % comp: Component class ref
        function comp = setVolMassTransFuncs(comp)
            % bub to bulk mass transfer
            % see line 1770 in ODESys.m
            % ### UPDATE: need better symbol for gas phase mol fraction
            % ### UPDATE: only add this function if necessary (if there is
            % an input of a gas as a bubble)
            % ### IMPT: this is done in ODESys.updateIOFunc()
            % bub_to_bulk = "(R*T)/((V_m-V)*"+comp.MW_sym+")*(D_i-(P*Y_i_"+comp.sym+"/"+comp.h_const_sym+"-"+comp.sym+"_bubf)/tau*V)"; % should be a mixing equation, not a transfer equation
            % 
            % % bub to liq mass transfer - ### FIXME: this needs improvement
            % % to better model how a PFR would work
            % bub_to_liq = "(P*Y_i_"+comp.sym+"/"+comp.h_const_sym+"-"+comp.sym+"_bub0)/tau_"+comp.sym; % discrete transfer rates calculated analytically
            % 
            % % bulk to liq mass transfer - based on concentration gradient
            % % at liquid surface/bulk gas interface
            % bulk_to_liq = "k_"+comp.sym+"_bulk*("+comp.bulk_gas_sym+"/"+comp.h_const_sym+"-"+comp.sym+")";
            % liq_to_bulk = "k_"+comp.sym+"_bulk*"+comp.h_const_sym+"*("+comp.sym+"-"+comp.bulk_gas_sym+")";
            % 
            % % bulk gas phase funcParams
            % comp.gasBulkFuncParams{1} = struct( ...
            %     'funcVal',bub_to_bulk, ...
            %     'funcName',[char(comp.sym),' bub_to_bulk'], ...
            %     'params',struct([]) ...
            % );
            % comp.gasBulkFuncParams{2} = struct( ...
            %     'funcVal',liq_to_bulk, ...
            %     'funcName',[char(comp.sym),' liq_to_bulk'], ...
            %     'params',struct([]) ...
            % );
            % % liquid phase funcParams
            % comp.funcParams{2} = struct( ...
            %     'funcVal',bub_to_liq, ...
            %     'funcName',[char(comp.sym),' bub_to_liq'], ...
            %     'params',struct([]) ...
            % );
            % comp.funcParams{3} = struct( ...
            %     'funcVal',bulk_to_liq, ...
            %     'funcName',[char(comp.sym),' bulk_to_liq'], ...
            %     'params',struct([]) ...
            % );
            % 
            % % setting parameters
            % funcVals = [bub_to_bulk,liq_to_bulk,bub_to_liq,bulk_to_liq];
            % depVars = [depVars(:,2);{['C_',char(string(comp.number))]}];
            % for k=1:1:length(funcVals)
            %     symVarStrArr = string(findVars(char(funcVals(k))));
            %     varStr = string(intersect(symVarStrArr,depVars));
            %     paramStr = convertStringsToChars(string(setxor(symVarStrArr,varStr)));
            %     for l=1:1:length(paramStr)
            %         switch k
            %             case 1
            %                 %(locNum, name, sym, val, unit, funcObj)
            %                 comp.gasBulkFuncParams{1} = comp.createNewParam(l,paramStr{l},paramStr{l},1,'',comp.gasBulkFuncParams{1},true);
            %             case 2
            %                 comp.gasBulkFuncParams{2} = comp.createNewParam(l,paramStr{l},paramStr{l},1,'',comp.gasBulkFuncParams{2},true);
            %             case 3
            %                 comp.funcParams{1} = comp.createNewParam(l,paramStr{l},paramStr{l},1,'',comp.funcParams{1},true);
            %             case 4
            %                 comp.funcParams{2} = comp.createNewParam(l,paramStr{l},paramStr{l},1,'',comp.funcParams{2},true);
            %         end
            %     end
            % end
        end

        % comp: Comp class ref
        function HConstFunc = getHConstHelperFuncs(comp,varargin)
            if ~isempty(varargin)
                if strcmp(varargin{1},"h_const")
                    HConstFunc = comp.h_const_helper_funcs{1};
                elseif strcmp(varargin{1},"eq_with_liq")
                    HConstFunc = comp.h_const_helper_funcs{2};
                elseif strcmp(varargin{1},"eq_with_gas")
                    HConstFunc = comp.h_const_helper_funcs{3};
                end
            else
                HConstFunc = comp.h_const_helper_funcs;
            end
        end

        % comp: Comp class ref
        function massTransFuncs = getMassTransHelpers(comp)
            massTransFuncs = comp.gasBulkHelpers;
        end

        % comp: Comp class ref
        function phases = getSorpPhases(comp)
            phases = {};
            for k=1:1:length(comp.sorpFuncParams)
                if ~any(strcmp(phases,comp.sorpFuncParams{k}.solventSym))
                    phases{end+1} = comp.sorpFuncParams{k}.solventSym;
                end
            end
        end

        % comp: Comp class ref, is_vol: Boolean
        function setVol(comp,is_vol,h_const,h_const_u,dh_const,dh_const_u,defaultParamVals)
            % ### FIXME: need to create helper for gas phase partial
            % pressure (convert from concentration in the gas phase with
            % Peng-Robinson or Ideal gas?)
            comp.is_vol = is_vol;

            % perform operations to set or reset volatility
            % properties
            if is_vol
                comp.bulk_gas_sym = char(replace(comp.sym,"C","P"));

                comp.h_const_u = char(string(h_const_u));
                switch comp.h_const_u
                    case 'kPa*L/g'
                    case 'kPa*L/mol'
                        comp.h_const = h_const./comp.MW;
                    case 'kPa*m^3/kg'
                        comp.h_const = h_const;
                end
                comp.dh_const_u = char(string(dh_const_u));
                switch comp.dh_const_u
                    case 'K'
                        comp.dh_const = dh_const;
                    case 'C'
                        comp.dh_const = dh_const;
                    case 'F'
                        comp.dh_const = dh_const .* 5 ./ 9;
                end
                comp.h_const_sym = ['H_0_',char(comp.sym)];
                comp.dh_const_sym = ['A_',char(comp.sym)];
                comp.h_const_helper_funcs = {};
                comp.setHConstHelperFuncs(defaultParamVals);
            else
                comp.bulk_gas_sym = '';
                comp.h_const_u = '';
                comp.h_const = 1;
                comp.dh_const_u = '';
                comp.dh_const = 1;
                comp.h_const_sym = '';
                comp.dh_const_sym = '';
                comp.h_const_helper_funcs = {};
            end
        end

        % comp: comps class ref, prop: string, val: number | string
        function comp = setInitConc(comp,conc,unit)
            comp.initConc = conc;
            comp.initConcUnit = unit;
        end

        % comp: Component class ref
        function initConc = getInitConc(comp)
            initConc = comp.initConc;
        end

        % comp: Component class ref
        function initConcUnit = getInitConcUnit(comp)
            initConcUnit = comp.initConcUnit;
        end

        % comp: Component class ref, initConc: num, is_vol: boolean,
        % h_const: num, h_const_u: string, dh_const: num, dh_const_u:
        % string
        function updateComp(comp,initConc,initConcUnit,is_vol,MW,h_const,h_const_u,dh_const,dh_const_u,defaultParamVals)
            comp.setInitConc(initConc,initConcUnit);
            comp.MW = MW;
            comp.setVol(is_vol,h_const,h_const_u,dh_const,dh_const_u,defaultParamVals);
        end

        % comp: Component class ref
        function type = getType(comp)
            type = comp.type;
        end

        % comp: Component class ref
        function clearFuncParams(comp)
            comp.funcParams = {};
        end

        % comp: comps class ref, funcVal: string, funcName: string,
        % paramStr: string{}, varNames: {}, defaultParamVals: struct
        function comp = setModel(comp, funcVal, funcName, paramStr, defaultParamVals)
            if isempty(funcName), funcName = comp.funcParams{1}.funcName; end
            comp.funcParams{1} = comp.reviseFuncParamObj(funcVal, funcName, comp.funcParams{1});
            prevParams = {};
            for k=1:1:length(comp.funcParams{1}.params)
                prevParams{k,1} = comp.funcParams{1}.params{k}.sym;
                prevParams{k,2} = comp.funcParams{1}.params{k}.val;
            end
            comp.funcParams{1} = comp.removeParams(comp.funcParams{1});
            % sets parameter syms
            locNum = 1;
            offset = 0;
            for k=1:1:length(paramStr)
                if ~isempty(prevParams) && any(strcmp(prevParams(:,1),paramStr{k}))
                    paramVal = prevParams{k-offset,2};
                else
                    paramVal = 1;
                    offset = offset + 1;
                end
                comp.funcParams{1} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, paramVal, '', comp.funcParams{1}, true, defaultParamVals);
                locNum = locNum + 1;
            end
        end

        % comp: Component class ref, funcVal: string, funcName: string,
        % paramStr: string[], phase: string, defaultParamVals: struct
        function comp = addModel(comp,funcVal,funcName,paramStr,phase,defaultParamVals)
            if strcmp(phase,'Main'), phase = 'Liquid'; end
            % sets parameter syms
            locNum = 1;
            switch phase
                case 'Liquid'
                    comp.funcParams{end+1} = comp.createFuncParamObj(funcVal,funcName);
                    for k=1:1:length(paramStr)
                        comp.funcParams{end} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, 1, '', comp.funcParams{end}, true, defaultParamVals);
                        locNum = locNum + 1;
                    end

                case 'Gas'
                    comp.gasBulkFuncParams{end+1} = comp.createFuncParamObj(funcVal,funcName);
                    for k=1:1:length(paramStr)
                        comp.gasBulkFuncParams{end} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, 1, '', comp.gasBulkFuncParams{end}, true, defaultParamVals);
                        locNum = locNum + 1;
                    end

                case 'Suspended Solid'
                    comp.sorpFuncParams{end+1} = comp.createFuncParamObj(funcVal,funcName);
                    for k=1:1:length(paramStr)
                        comp.sorpFuncParams{end} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, 1, '', comp.sorpFuncParams{end}, true, defaultParamVals);
                        locNum = locNum + 1;
                    end
            end
        end

        % comp: Component class ref
        function comp = rmModel(comp)
            % input string of which funcParam group to remove from
            comp.removeFuncParamObj();
        end

        % comp: Component class ref, compName: string, funcVal: string,
        % funcName: string, paramStr: string{}, defaultParamVals: struct
        function comp = addSuspendedSolidVolumeFunc(comp,funcVal,funcName,paramStr,defaultParamVals)
            comp.sorpVolFuncParams{end+1} = comp.createFuncParamObj(funcVal,funcName);

            % sets parameter syms
            locNum = 1;
            for k=1:1:length(paramStr)
                comp.sorpVolFuncParams{end} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, 1, '', comp.sorpVolFuncParams{end}, true, defaultParamVals);
                locNum = locNum + 1;
            end
        end

        % comp: Component class ref, compName: string, funcVal: string,
        % funcName: string
        function comp = rmSuspendedSolidVolumeFunc(comp,funcVal,funcName)
            for k=1:1:length(comp.sorpVolFuncParams)
                if strcmp(comp.sorpVolFuncParams{k}.funcVal,funcVal) && ...
                        strcmp(comp.sorpVolFuncParams{k}.funcName,funcName)
                    comp.sorpVolFuncParams(k) = [];
                    return;
                end
            end
        end

        % comp: Component class ref, funcVal: string, phaseA: string,
        % phaseB: string, paramStr: string{}, defaultParamVals: struct
        function comp = setMTModel(comp,funcVal,phaseA,phaseB,funcType,paramStr,defaultParamVals)
            % updating the function values in the correct funcParam
            if strcmp(phaseA,'Liquid')
                idx = 0;
                for k=1:1:length(comp.funcParams)
                    if strcmp(comp.funcParams{k}.funcName,[phaseB,' MT'])
                        idx = k;
                    end
                end

                prevParams = {};
                for k=1:1:length(comp.funcParams{idx}.params)
                    prevParams{k,1} = comp.funcParams{idx}.params{k}.sym;
                    prevParams{k,2} = comp.funcParams{idx}.params{k}.val;
                end

                comp.funcParams{idx} = comp.reviseFuncParamObj(funcVal,[phaseB,' MT'],comp.funcParams{idx},funcType);
                comp.funcParams{idx} = comp.removeParams(comp.funcParams{idx});
                % sets parameter syms
                locNum = 1;
                for k=1:1:length(paramStr)
                    if ~isempty(prevParams) && any(strcmp(prevParams(:,1),paramStr{k}))
                        paramVal = prevParams{k,2};
                    else
                        paramVal = 1;
                    end
                    comp.funcParams{idx} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, paramVal, '', comp.funcParams{idx}, true, defaultParamVals);
                    locNum = locNum + 1;
                end
            elseif strcmp(phaseA,'Gas')
                idx = 1;
                for k=1:1:length(comp.gasBulkFuncParams)
                    if strcmp(comp.gasBulkFuncParams{k}.funcName,[phaseB,' MT'])
                        idx = k;
                    end
                end
                
                prevParams = {};
                for k=1:1:length(comp.gasBulkFuncParams{idx}.params)
                    prevParams{k,1} = comp.gasBulkFuncParams{idx}.params{k}.sym;
                    prevParams{k,2} = comp.gasBulkFuncParams{idx}.params{k}.val;
                end

                comp.gasBulkFuncParams{idx} = comp.reviseFuncParamObj(funcVal,[phaseB,' MT'],comp.gasBulkFuncParams{idx},funcType);
                comp.gasBulkFuncParams{idx} = comp.removeParams(comp.gasBulkFuncParams{idx});
                % sets parameter syms
                locNum = 1;
                for k=1:1:length(paramStr)
                    if ~isempty(prevParams) && any(strcmp(prevParams(:,1),paramStr{k}))
                        paramVal = prevParams{k,2};
                    else
                        paramVal = 1;
                    end
                    comp.gasBulkFuncParams{idx} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, paramVal, '', comp.gasBulkFuncParams{end}, true, defaultParamVals);
                    locNum = locNum + 1;
                end
            else
                idx = 1;
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseA)
                        idx = k;
                    end
                end
                comp.sorpFuncParams{idx} = comp.reviseFuncParamObj(funcVal,phaseA,comp.sorpFuncParams{idx},funcType);
                locNum = 1;
                for k=1:1:length(paramStr)
                    comp.sorpFuncParams{idx} = comp.createNewParam(locNum, paramStr{k}, paramStr{k}, 1, '', comp.sorpFuncParams{end}, true, defaultParamVals);
                    locNum = locNum + 1;
                end
            end
        end

        % comp: comps class ref, chemName: string, func: string
        function funcExists = checkForDuplicateFunc(comp, funcVal, funcName)
            funcExists = false;
            for k=1:1:numel(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcVal,funcVal) && ...
                        strcmp(comp.funcParams{k}.funcName,funcName)
                    funcExists = true;
                    return;
                end
            end
        end

        % comp: comps class ref
        function sym = getSym(comp)
            sym = comp.sym;
        end

        % comp: Component class ref
        function sym = getBulkGasSym(comp)
            sym = char(comp.bulk_gas_sym);
        end

        % comp: Component class ref
        function sym = getMWSym(comp)
            sym = char(comp.MW_sym);
        end

        % comp: Component class ref
        function name = getName(comp)
            name = comp.name;
        end

        % comp: Component class ref
        function names = getAllPhaseNames(comp)
            names = comp.getName()+" (Liquid Phase)";
            if comp.is_vol
                names(end+1) = comp.getName()+" (Gas Phase)";
            end
            for k=1:1:length(comp.sorpFuncParams)
                names(end+1) = comp.getName()+" ("+comp.sorpFuncParams{k}.solventName+" Phase)";
            end
        end

        % comp: Component class ref
        function num = getNum(comp)
            num = comp.number;
        end

        % comp: comps class ref, phase: string, phaseName: string, dispMT: boolean
        function funcs = getGovFuncs(comp, phase, phaseName)
            funcs = {};
            if strcmp(phase,'Gas Phase')
                fp = comp.gasBulkFuncParams;
            elseif strcmp(phase,'Liquid Phase')
                fp = comp.funcParams;
            elseif strcmp(phase,'Suspended Solid Phase')
                fp = {};
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseName)
                        fp{end+1} = comp.sorpFuncParams{k}; %#ok<AGROW>
                    end
                end
            end
            for k=1:1:length(fp)
                % converting piecewise expressions to latex
                funcVal = comp.convertPWToLaTeX(fp{k});
                funcs{end+1} = funcVal; %#ok<AGROW>
            end
            if isempty(funcs)
                funcs{1} = "0";
            end
        end

        % comp: comps class ref, field: string[], val: cell{}
        % returns the indices for comp.funcParams that correspond to items
        % that match the fields searched
        % make sure that val inputs for text are strings
%         function idx = getFuncByField(comp, fields, val)
%             funcObjs1 = {};
%             for k=1:1:length(comp.funcParams)
%                 for l=1:1:length(fieldnames(comp.funcParams{1}))
%                     funcObjs2 = {};
%                     if comp.funcParams{k}.(field{l}) == val{l}
%                         funcObjs2{l} = comp.funcParams{k}; %#ok<AGROW> 
%                     end
%                     funcObjs1{k} = funcObjs2; %#ok<AGROW> 
%                 end
%             end
% 
%             funcObjs = intersect(funcObjs1{1},funcObjs1{2});
%             for k=3:1:length(funcObjs1)
%                 funcObjs = intersect(funcObjs,funcObjs1{k});
%             end

%             comboFuncParam = comp.funcParams{:};
%             idxArr = zeros(length(fields),length(comp.funcParams))
%             disp(comboFuncParam)
%             disp([comboFuncParam.(fields(1))])
%             disp(val{1})
%             for k=1:1:length(fields)
%                 if isa(val{k},'char')
%                     val{k} = convertCharsToStrings(val{k});
%                 end
%                 disp([comboFuncParam.(fields(k))])
%                 idxArr(k,:) = ([comboFuncParam.(fields(k))] == val{k});
%             end
%             if length(fields) < 2
%                 idx = idxArr;
%             end
%             for k=1:1:(length(fields)-1)
%                 idx = intersect(idxArr(k,:),idxArr(k+1,:));
%             end
%         end

        % comp: comps class ref, paramSym: string, newVal: number,
        % newUnit: string, newParamName: string, funcName: string,
        % editable: boolean
        function comp = updateParam(comp, paramSym, newVal, newUnit, newParamName, funcName, convertUnits)
            try
                funcObjIdx = comp.getFuncIdx(funcName);
            catch
                disp("err Component.m Line 412")
                return;
            end
            if convertUnits
                newVal = unit_standardization(newVal,newUnit);
            end

            comp.funcParams{funcObjIdx} = comp.reviseParam(paramSym,newVal,newUnit,newParamName,comp.funcParams{funcObjIdx});
        end

        % comp: comps class ref, funcName: string,
        function idx = getFuncIdx(comp, funcName)
            % ### FIXME: if functions have the same value, then there is an
            % error because the first matching function is pulled
            for k=1:1:length(comp.funcParams)
                if strcmp(comp.funcParams{k}.funcName,funcName)
                    idx = k;
                    return;
                end
            end
            disp("No Function Found, Component.m, line 223");
        end

        % comp: comps class ref
        % returns parameter info ready to display to ODESys (params struct)
        function params = getGrthParams(comp,phase)
            if strcmp(phase,'Main'), phase = 'Liquid'; end
            params = cell([comp.getParamCount(phase),6]);
            ct = 1;
            if strcmp(phase,'Liquid')
                fp = comp.funcParams;
            elseif strcmp(phase,'Gas')
                fp = comp.gasBulkFuncParams;
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(phase,comp.sorpFuncParams{k}.solventName)
                        fp = comp.sorpFuncParams{k};
                    end
                end
            end
            for k=1:1:length(fp)
                for l=1:1:length(fp{k}.params)
                    val = unit_standardization_inv(fp{k}.params{l}.val,fp{k}.params{l}.unit);

                    params{ct,1} = num2str(ct);
                    params{ct,2} = fp{k}.funcName;
                    params{ct,3} = fp{k}.params{l}.sym;
                    params{ct,4} = fp{k}.params{l}.name;
                    params{ct,5} = val;
                    params{ct,6} = fp{k}.params{l}.unit;
                    params{ct,7} = fp{k}.params{l}.editable;
                    ct = ct + 1;
                end
            end
        end

        % comp: comps class ref, phase: string
        function paramCt = getParamCount(comp,phase)
            paramCt = 0;
            if strcmp(phase,'Liquid')
                fp = comp.funcParams;
            elseif strcmp(phase,'Gas')
                fp = comp.gasBulkFuncParams;
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(phase,comp.sorpFuncParams{k}.solventName)
                        fp = comp.sorpFuncParams{k};
                    end
                end
            end
            for k=1:1:length(fp)
                paramCt = paramCt + length(fp{k}.params);
            end
        end

        % comp: comps class ref
        function removeModel(comp,funcVal,funcName,phase)
            if strcmp(phase,'Main'), phase = 'Liquid'; end
            switch phase
                case 'Liquid'
                    for k=1:1:length(comp.funcParams)
                        if strcmp(comp.funcParams{k}.funcVal,funcVal) && ...
                                strcmp(comp.funcParams{k}.funcName,funcName)
                            comp.funcParams(k) = [];
                            return;
                        end
                    end
                    if isempty(comp.funcParams)
                        
                    end

                case 'Gas'
                    for k=1:1:length(comp.gasBulkFuncParams)
                        if strcmp(comp.gasBulkFuncParams{k}.funcVal,funcVal) && ...
                                strcmp(comp.gasBulkFuncParams{k}.funcName,funcName)
                            comp.gasBulkFuncParams(k) = [];
                            return;
                        end
                    end

                case 'Suspended Solid'
                    for k=1:1:length(comp.sorpFuncParams)
                        if strcmp(comp.sorpFuncParams{k}.funcVal,funcVal) && ...
                                strcmp(comp.sorpFuncParams{k}.funcName,funcName)
                            comp.sorpFuncParams(k) = [];
                            return;
                        end
                    end
            end
        end

        % comp: Component class ref, funcName: string, phase: string
        function removeModelByName(comp,funcName,phase)
            if strcmp(phase,'Main'), phase = 'Liquid'; end
            switch phase
                case 'Liquid'
                    for k=1:1:length(comp.funcParams)
                        if strcmp(comp.funcParams{k}.funcName,funcName)
                            comp.funcParams(k) = [];
                            return;
                        end
                    end

                case 'Gas'
                    for k=1:1:length(comp.gasBulkFuncParams)
                        if strcmp(comp.gasBulkFuncParams{k}.funcName,funcName)
                            comp.gasBulkFuncParams(k) = [];
                            return;
                        end
                    end

                case 'Suspended Solid'
                    for k=1:1:length(comp.sorpFuncParams)
                        if strcmp(comp.sorpFuncParams{k}.funcName,funcName)
                            comp.sorpFuncParams(k) = [];
                            return;
                        end
                    end
            end
        end

        % comp: Component class ref
        function funcNum = getFuncCount(comp)
            funcNum = length(comp.funcParams);
        end

        % comp: Component class ref, reg_param_ct: number, paramCt: numbers, paramList: {}, regParamList: {}, regression: boolean
        function [govFunc,params,reg_param_ct,regParamListUpdate] = compileGovFunc(comp,param_ct,reg_param_ct,regParamList,regression)
            % Need to make sure that each individual govFunc component is
            % isolated (so (C1+C2)*(C3+C4) ~= C1+C2*C3+C4)
            regParamListUpdate = {};
            govFunc = "";
            params = [];

            if strcmp(comp.type,'Suspended Solid Sorbent')
                % compiling governing functions (biological and mass transfer)
                for k=1:1:length(comp.sorpVolFuncParams)
                    % compilation for piecewise functions
                    funcVal = comp.compilePWFuncs(comp.sorpVolFuncParams{k});
                    charGovFunc = char(govFunc);
                    if isempty(charGovFunc) || strcmp(charGovFunc(end),';')
                        operator = "";
                    else
                        operator = "+";
                    end
                    govFunc = govFunc + operator + funcVal;
                    % compiling parameters array
                    for l=1:1:length(comp.sorpVolFuncParams{k}.params)
                        params = [params,comp.sorpVolFuncParams{k}.params{l}.val];
                    end
                end
                fp = comp.sorpVolFuncParams;
            else
                % compiling governing functions (biological and mass transfer)
                for k=1:1:length(comp.funcParams)
                    % compilation for piecewise functions
                    funcVal = comp.compilePWFuncs(comp.funcParams{k});
                    charGovFunc = char(govFunc);
                    if isempty(charGovFunc) || strcmp(charGovFunc(end),';')
                        operator = "";
                    else
                        operator = "+";
                    end
                    govFunc = govFunc + operator + funcVal;
                    % compiling parameters array
                    for l=1:1:length(comp.funcParams{k}.params)
                        params = [params,comp.funcParams{k}.params{l}.val];
                    end
                end
                % Gas concentration - should only need to track the bulk
                % concentration, the bubble concentration will be an analytical
                % calculation, add the mass transfer functions from perspective
                % of bulk, bubble mass addition and mass transfer between bulk
                % and liquid
                if comp.is_vol
                    govFunc = govFunc + ";";
                    % compilation for piecewise functions
                    % compiling governing functions (biological and mass transfer)
                    for k=1:1:length(comp.gasBulkFuncParams)
                        funcVal = comp.compilePWFuncs(comp.gasBulkFuncParams{k});
                        charGovFunc = char(govFunc);
                        if isempty(charGovFunc) || strcmp(charGovFunc(end),';')
                            operator = "";
                        else
                            operator = "+";
                        end
                        govFunc = govFunc + operator + funcVal;
                        % compiling parameters array
                        for l=1:1:length(comp.gasBulkFuncParams{k}.params)
                            params = [params,comp.gasBulkFuncParams{k}.params{l}.val];
                        end
                    end
                end
                % compiling sorption phases
                if ~isempty(comp.sorpFuncParams)
                    govFunc = govFunc + ";";
                    % compilation for piecewise functions
                    for k=1:1:length(comp.sorpFuncParams)
                        funcVal = comp.compilePWFuncs(comp.sorpFuncParams{k});
                        charGovFunc = char(govFunc);
                        if isempty(charGovFunc) || strcmp(charGovFunc(end),';')
                            operator = "";
                        else
                            operator = "+";
                        end
                        govFunc = govFunc + operator + funcVal;
                        % compiling parameters array
                        for l=1:1:length(comp.sorpFuncParams{k}.params)
                            params = [params,comp.sorpFuncParams{k}.params{l}.val];
                        end
                    end
                end
                fp = [comp.funcParams,comp.gasBulkFuncParams,comp.sorpFuncParams];
            end

            % replacing param syms in govFunc with p({gloNum})
            ct = 0;
            reg_param_ct_loc = 0;
            if regression
                for k=1:1:length(fp)
                    % replacing function and regression parameters that haven't yet been
                    % replaced
                    for l=1:1:length(fp{k}.params)
                        match = strcmp(regParamList(:,1),char(fp{k}.params{l}.sym));
                        % ### NOTE: replace variables passed in
                        % in Simulink as v (like p but v)
                        ct = ct + 1;
                        if any(match)
                            if ~strcmp(regParamList{match,2},"")
                                govFunc = replace(govFunc,fp{k}.params{l}.sym,"$("+(regParamList{match,2})+")");
                            else
                                reg_param_ct_loc = reg_param_ct_loc + 1;
                                govFunc = replace(govFunc,fp{k}.params{l}.sym,"$("+(reg_param_ct+reg_param_ct_loc)+")");
                                regParamListUpdate{1,match} = true; %#ok<*AGROW>
                                regParamListUpdate{2,match} = find(match);
                                regParamListUpdate{3,match} = string(reg_param_ct+reg_param_ct_loc);
                            end
                        else
                            govFunc = replace(govFunc,fp{k}.params{l}.sym,"#("+(param_ct+ct)+")");
                        end
                    end
                end
            else
                % replacing function parameters that haven't been replaced
                % yet
                for k=1:1:length(fp)
                    for l=1:1:length(fp{k}.params)
                        ct = ct + 1;
                        govFunc = replace(govFunc,fp{k}.params{l}.sym,"#("+(param_ct+ct)+")");
                    end
                end
            end
            % govFunc = strrep(strrep(strjoin(govFunc),"#","")," ","");
            reg_param_ct = reg_param_ct + reg_param_ct_loc;
            govFunc = replace(replace(govFunc,"#","p"),"$","r");
            if strcmp(govFunc,""), govFunc = "0"; end
        end

        % comp: Compnent class ref, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, lowLimOp: string,
        % upLimOp: string, phaseA: string, phaseB: string
        function pwTInfo = addNewPWExpr(comp,newExpr,lowLim,upLim,sysVar,elseVal,lowLimOp,upLimOp,phaseA,phaseB)
            if strcmp(phaseA,'Main')
                comp.funcParams{1}.lims{end+1} = struct([ ...
                    'funcVal',newExpr, ...
                    'lowerLim',lowLim ...
                    'upperLim',upLim, ...
                    'sysVar',sysVar, ...
                    'elseVal',elseVal, ...
                    'lowLimOp',lowLimOp, ...
                    'upLimOp',upLimOp
                ]);
                pwTInfo = cell(length(comp.funcParams{1}.lims),4);
                for k=1:1:length(comp.funcParams{1}.lims)
                    pwTInfo{k,1} = comp.funcParams{1}.lims{k}.lowerLim;
                    pwTInfo{k,2} = comp.funcParams{1}.lims{k}.upperLim;
                    pwTInfo{k,3} = comp.funcParams{1}.lims{k}.sysVar;
                    pwTInfo{k,4} = comp.funcParams{1}.lims{k}.elseVal;
                end
            elseif strcmp(phaseA,'Liquid')
                for k=2:1:length(comp.funcParams)
                    if strcmp(comp.funcParams{k}.funcName,[phaseB,' MT'])
                        comp.funcParams{k}.lims{end+1} = struct([ ...
                            'funcVal',newExpr, ...
                            'lowerLim',lowLim ...
                            'upperLim',upLim, ...
                            'sysVar',sysVar, ...
                            'elseVal',elseVal, ...
                            'lowLimOp',lowLimOp, ...
                            'upLimOp',upLimOp
                        ]);
                        pwTInfo = cell(length(comp.funcParams{k}.lims),4);
                        for l=1:1:length(comp.funcParams{k}.lims)
                            pwTInfo{l,1} = comp.funcParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.funcParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.funcParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.funcParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            elseif strcmp(phaseA,'Gas')
                for k=1:1:length(comp.gasBulkFuncParams)
                    if strcmp(comp.gasBulkFuncParams{k}.funcName,[phaseB,' MT'])
                        comp.gasBulkFuncParams{k}.lims{end+1} = struct([ ...
                            'funcVal',newExpr, ...
                            'lowerLim',lowLim ...
                            'upperLim',upLim, ...
                            'sysVar',sysVar, ...
                            'elseVal',elseVal, ...
                            'lowLimOp',lowLimOp, ...
                            'upLimOp',upLimOp
                        ]);
                        pwTInfo = cell(length(comp.gasBulkFuncParams{k}.lims),4);
                        for l=1:1:length(comp.gasBulkFuncParams{k}.lims)
                            pwTInfo{l,1} = comp.gasBulkFuncParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.gasBulkFuncParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.gasBulkFuncParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.gasBulkFuncParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseA)
                        comp.sorpFuncParams{k}.lims{end+1} = struct([ ...
                            'funcVal',newExpr, ...
                            'lowerLim',lowLim ...
                            'upperLim',upLim, ...
                            'sysVar',sysVar, ...
                            'elseVal',elseVal, ...
                            'lowLimOp',lowLimOp, ...
                            'upLimOp',upLimOp
                        ]);
                        pwTInfo = cell(length(comp.sorpFuncParams{k}.lims),4);
                        for l=1:1:length(comp.sorpFuncParams{k}.lims)
                            pwTInfo{l,1} = comp.sorpFuncParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.sorpFuncParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.sorpFuncParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.sorpFuncParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            end
        end

        % comp: Compnent class ref, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, lowLimOp: string,
        % upLimOp: string, phaseA: string, phaseB: string
        function pwTInfo = editPWExpr(comp,newExpr,lowLim,upLim,sysVar,elseVal,phaseA,phaseB)
            idx = 1;
            newLim = struct([ ...
                    'funcVal',newExpr, ...
                    'lowerLim',lowLim ...
                    'upperLim',upLim, ...
                    'sysVar',sysVar, ...
                    'elseVal',elseVal, ...
                    'lowLimOp',lowLimOp, ...
                    'upLimOp',upLimOp
            ]);
            if strcmp(phaseA,'Main')
                for k=1:1:length(comp.funcParams{1}.lims)
                    if strcmp(comp.funcParams{1}.lims{k}.funcVal,newExpr)
                        idx = k;
                        break;
                    end
                end
                comp.funcParams{1}.lims{idx} = newLim;
                pwTInfo = cell(length(comp.funcParams{1}.lims),4);
                for k=1:1:length(comp.funcParams{1}.lims)
                    pwTInfo{k,1} = comp.funcParams{1}.lims{k}.lowerLim;
                    pwTInfo{k,2} = comp.funcParams{1}.lims{k}.upperLim;
                    pwTInfo{k,3} = comp.funcParams{1}.lims{k}.sysVar;
                    pwTInfo{k,4} = comp.funcParams{1}.lims{k}.elseVal;
                end
            elseif strcmp(phaseA,'Liquid')
                for k=1:1:length(comp.funcParams)
                    if strcmp(comp.funcParams{k}.solventName,phaseB)
                        for l=1:1:length(comp.funcParams{k}.lims)
                            if strcmp(comp.funcParams{k}.lims{l}.funcVal,newExpr)
                                idx = l;
                                break;
                            end
                        end
                        comp.funcParams{k}.lims{idx} = newLim;
                        pwTInfo = cell(length(comp.funcParams{1}.lims),4);
                        for l=1:1:length(comp.funcParams{k}.lims)
                            pwTInfo{l,1} = comp.funcParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.funcParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.funcParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.funcParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            elseif strcmp(phaseA,'Gas')
                for k=1:1:length(comp.gasBulkFuncParams)
                    if strcmp(comp.gasBulkFuncParams{k}.funcName,[phaseB,' MT'])
                        for l=1:1:length(comp.gasBulkFuncParams{k}.lims)
                            if strcmp(comp.gasBulkFuncParams{k}.lims{l}.funcVal,newExpr)
                                idx = l;
                                break;
                            end
                        end
                        comp.gasBulkFuncParams{k}.lims{idx} = newLim;
                        pwTInfo = cell(length(comp.gasBulkFuncParams{1}.lims),4);
                        for l=1:1:length(comp.gasBulkFuncParams{k}.lims)
                            pwTInfo{l,1} = comp.gasBulkFuncParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.gasBulkFuncParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.gasBulkFuncParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.gasBulkFuncParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseA)
                        for l=1:1:length(comp.gasBulkFuncParams{k}.lims)
                            if strcmp(comp.sorpFuncParams{k}.lims{l}.funcVal,newExpr)
                                idx = l;
                                break;
                            end
                        end
                        comp.sorpFuncParams{k}.lims{idx} = newLim;
                        pwTInfo = cell(length(comp.sorpFuncParams{k}.lims),4);
                        for l=1:1:length(comp.sorpFuncParams{k}.lims)
                            pwTInfo{l,1} = comp.sorpFuncParams{k}.lims{l}.lowerLim;
                            pwTInfo{l,2} = comp.sorpFuncParams{k}.lims{l}.upperLim;
                            pwTInfo{l,3} = comp.sorpFuncParams{k}.lims{l}.sysVar;
                            pwTInfo{l,4} = comp.sorpFuncParams{k}.lims{l}.elseVal;
                        end
                    end
                end
            end
        end

        % comp: Compnent class ref, newExpr: string, lowLim: number, 
        % upLim: number, sysVar: string, elseVal: number, phaseA: string,
        % phaseB: string
        function pwTInfo = rmPWExpr(comp,newExpr,lowLim,upLim,sysVar,elseVal,phaseA,phaseB)
            idx = 1;
            if strcmp(phaseA,'Main')
                for k=1:1:length(comp.funcParams{1}.lims)
                    if strcmp(comp.funcParams{1}.lims{k}.funcVal,newExpr) && ...
                            comp.funcParams{1}.lims{k}.lowerLim == lowLim && ...
                            comp.funcParams{1}.lims{k}.upperLim == upLim && ...
                            strcmp(comp.funcParams{1}.lims{k}.sysVar,sysVar) && ...
                            strcmp(comp.funcParams{1}.lims{k}.elseVal,elseVal)
                        idx = k;
                        break;
                    end
                end
                comp.funcParams{1}.lims(idx) = [];
                pwTInfo = cell(length(comp.funcParams{1}.lims),4);
                for k=1:1:length(comp.funcParams{1}.lims)
                    pwTInfo{k,1} = comp.funcParams{1}.lims{k}.lowerLim;
                    pwTInfo{k,2} = comp.funcParams{1}.lims{k}.upperLim;
                    pwTInfo{k,3} = comp.funcParams{1}.lims{k}.sysVar;
                    pwTInfo{k,4} = comp.funcParams{1}.lims{k}.elseVal;
                end
            elseif strcmp(phaseA,'Liquid')
                for k=1:1:length(comp.funcParams)
                    if strcmp(comp.funcParams{k}.solventName,phaseB)
                        for l=1:1:length(comp.funcParams{k}.lims)
                            if strcmp(comp.funcParams{k}.lims{l}.funcVal,newExpr) && ...
                                    comp.funcParams{k}.lims{l}.lowerLim == lowLim && ...
                                    comp.funcParams{k}.lims{l}.upperLim == upLim && ...
                                    strcmp(comp.funcParams{k}.lims{l}.sysVar,sysVar) && ...
                                    strcmp(comp.funcParams{k}.lims{l}.elseVal,elseVal)
                                idx = l;
                                break;
                            end
                        end
                    end
                end
                comp.funcParams{k}.lims(idx) = [];
                pwTInfo = cell(length(comp.funcParams{k}.lims),4);
                for k=1:1:length(comp.funcParams{k}.lims)
                    pwTInfo{l,1} = comp.funcParams{k}.lims{l}.lowerLim;
                    pwTInfo{l,2} = comp.funcParams{k}.lims{l}.upperLim;
                    pwTInfo{l,3} = comp.funcParams{k}.lims{l}.sysVar;
                    pwTInfo{l,4} = comp.funcParams{k}.lims{l}.elseVal;
                end
            elseif strcmp(phaseA,'Gas')
                for k=1:1:length(comp.gasBulkFuncParams)
                    if strcmp(comp.gasBulkFuncParams{k}.funcName,[phaseB,' MT'])
                        for l=1:1:length(comp.gasBulkFuncParams{k}.lims)
                            if strcmp(comp.gasBulkFuncParams{k}.lims{l}.funcVal,newExpr) && ...
                                    comp.gasBulkFuncParams{k}.lims{l}.lowerLim == lowLim && ...
                                    comp.gasBulkFuncParams{k}.lims{l}.upperLim == upLim && ...
                                    strcmp(comp.gasBulkFuncParams{k}.lims{l}.sysVar,sysVar) && ...
                                    strcmp(comp.gasBulkFuncParams{k}.lims{l}.elseVal,elseVal)
                                idx = l;
                                break;
                            end
                        end
                    end
                end
                comp.gasBulkFuncParams{k}.lims(idx) = [];
                pwTInfo = cell(length(comp.gasBulkFuncParams{k}.lims),4);
                for k=1:1:length(comp.gasBulkFuncParams{k}.lims)
                    pwTInfo{l,1} = comp.gasBulkFuncParams{k}.lims{l}.lowerLim;
                    pwTInfo{l,2} = comp.gasBulkFuncParams{k}.lims{l}.upperLim;
                    pwTInfo{l,3} = comp.gasBulkFuncParams{k}.lims{l}.sysVar;
                    pwTInfo{l,4} = comp.gasBulkFuncParams{k}.lims{l}.elseVal;
                end
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseA)
                        for l=1:1:length(comp.sorpFuncParams{k}.lims)
                            if strcmp(comp.sorpFuncParams{k}.lims{l}.funcVal,newExpr) && ...
                                    comp.sorpFuncParams{k}.lims{l}.lowerLim == lowLim && ...
                                    comp.sorpFuncParams{k}.lims{l}.upperLim == upLim && ...
                                    strcmp(comp.sorpFuncParams{k}.lims{l}.sysVar,sysVar) && ...
                                    strcmp(comp.sorpFuncParams{k}.lims{l}.elseVal,elseVal)
                                idx = l;
                                break;
                            end
                        end
                    end
                end
                comp.sorpFuncParams{k}.lims(idx) = [];
                pwTInfo = cell(length(comp.sorpFuncParams{k}.lims),4);
                for k=1:1:length(comp.sorpFuncParams{k}.lims)
                    pwTInfo{l,1} = comp.sorpFuncParams{k}.lims{l}.lowerLim;
                    pwTInfo{l,2} = comp.sorpFuncParams{k}.lims{l}.upperLim;
                    pwTInfo{l,3} = comp.sorpFuncParams{k}.lims{l}.sysVar;
                    pwTInfo{l,4} = comp.sorpFuncParams{k}.lims{l}.elseVal;
                end
            end
        end

        % comp: Component class ref, solventName: string, solventNum: num,
        % soluteName: string, soluteNum: num,
        % modelName: string, param1Val: num, param2Val: num, param3Val: num
        function res = addSorptionEq(comp,solventName,solventNum,soluteName,soluteNum,modelName,param1Val,param2Val,param3Val,defaultParamVals)
            % add functionality to track info for sorbent
            % modeling
            if strcmp(comp.type,'Chemical Solute')
                funcVal = "0";
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,solventName)
                        funcVal = comp.sorpFuncParams{k}.funcVal;
                        comp.rmSorptionEq(solventName,solventNum,soluteName,modelName,param1Val,param2Val,param3Val);
                    end
                end
                comp.createSorpFuncParam(solventName,solventNum,modelName,funcVal,param1Val,param2Val,param3Val,defaultParamVals);
                res = comp.sorpHelpers;
            elseif strcmp(comp.type,'Suspended Solid Sorbent')
                if size(comp.sorpTInfo,2) > 0
                    for k=1:1:length(comp.sorpTInfo(:,1))
                        if strcmp(comp.sorpTInfo(:,1),soluteName)
                            comp.rmSorptionEq(solventName,solventNum,soluteName,modelName,param1Val,param2Val,param3Val);
                        end
                    end
                end
                comp.sorpTInfo{end+1,1} = soluteName;
                comp.sorpTInfo{end,2} = modelName;
                switch modelName
                    case 'Partition Coefficient'
                        comp.sorpTInfo{end,3} = ['K_C_',char(string(soluteNum)),'_',char(string(solventNum)),'=',char(string(param1Val))];

                    case 'Freundlich Adsorption'


                end

                res = comp.sorpTInfo;
            end
        end

        % comp: Component class ref, solventName: string, solventNum: num,
        % soluteName: string, soluteNum: num,
        % modelName: string, param1Val: num, param2Val: num, param3Val: num
        function res = rmSorptionEq(comp,solventName,~,soluteName,~,~,~,~)
            if strcmp(comp.type,'Chemical Solute')
                comp.removeSorpFuncParamObj(solventName);
                res = {};
                for k=1:1:length(comp.sorpHelpers)
                    if contains(comp.sorpHelpers{k}.funcName,solventName) 
                        res{end+1} = comp.sorpHelpers{k}; %#ok<AGROW>
                    end
                end
            elseif strcmp(comp.type,'Suspended Solid Sorbent')
                for k=1:1:size(comp.sorpTInfo,1)
                    if strcmp(comp.sorpTInfo{k,1},soluteName)
                        comp.sorpTInfo(k,:) = [];
                        res = comp.sorpTInfo;
                        return;
                    end
                end
            end
        end

        % comp: Component class ref, solventName: string, solventNum: num, funcVal: string,
        % modeName: string funcVal: string, param1Val: num, param2Val, numm, param3Val: num
        function comp = createSorpFuncParam(comp,solventName,solventNum,modelName,funcVal,param1Val,param2Val,param3Val,defaultParamVals)
            % creates funcParams for sorption MT governing
            % functions
            funcObj = struct( ...
                'solventName',solventName, ...
                'solventSym',['S_',char(string(comp.number)),'_',char(string(solventNum))], ...
                'funcName',[solventName,' MT'], ...
                'funcType','None', ...
                'funcVal',funcVal, ...
                'params',struct([]), ...
                'lims',cell(1,1) ...
            );
    
            switch modelName
                case 'Partition Coefficient'
                    % create helper functions for liquid phase concentration in
                    % equilibrium with sorped phase concentration, and
                    % sorped phase concentration in equilibrium with liquid
                    % phase concentration
                    partCoeffFuncVal = ['K_C0_',char(string(comp.number)),'_',char(string(solventNum))];
                    comp.sorpHelpers{end+1} = SubFunc(partCoeffFuncVal,[char(comp.name),' Partition Coefficient between Liquid Phase and ',solventName], ...
                        ['K_C_',char(string(comp.number)),'_',char(string(solventNum))],0,0);
                    comp.sorpHelpers{end}.updateParams([char(comp.name),' Partition Coefficient between Liquid Phase and ',solventName], ...
                        ['K_C0_',char(string(comp.number)),'_',char(string(solventNum))],param1Val,'-',false,defaultParamVals);

                    eqLiqConcFuncVal = ['S_',char(string(comp.number)),'_',char(string(solventNum)),'/K_C_',char(string(comp.number)),'_',char(string(solventNum))];
                    comp.sorpHelpers{end+1} = SubFunc(eqLiqConcFuncVal,[char(comp.name),' Concentration in Liquid Phase in Equilibrium with ',solventName], ...
                        ['C_eq_',char(string(comp.number)),'_',char(string(solventNum))],0,0);
                    comp.sorpHelpers{end}.updateParams([char(comp.name),' Partition Coefficient between Liquid Phase and ',solventName], ...
                        ['K_C_',char(string(comp.number)),'_',char(string(solventNum))],param1Val,'-',false,defaultParamVals);
                    
                    eqSorpConcFuncVal = ['K_C_',char(string(comp.number)),'_',char(string(solventNum)),'*C_',char(string(comp.number))];
                    comp.sorpHelpers{end+1} = SubFunc(eqSorpConcFuncVal,[char(comp.name),' Concentration in ',solventName,' Phase in Equilibrium with Liquid Phase'], ...
                        ['S_eq_',char(string(comp.number)),'_',char(string(solventNum))],0,0);
                    comp.sorpHelpers{end}.updateParams([char(comp.name),' Partition Coefficient between Liquid Phase and ',solventName], ...
                        ['K_C_',char(string(comp.number)),'_',char(string(solventNum))],param1Val,'-',false,defaultParamVals);

                case 'Freundlich Adsorption'
                    % create helper functions for liquid phase concentration in
                    % equilibrium with sorped phase concentration, and
                    % sorped phase concentration in equilibrium with liquid
                    % phase concentration


                    % creates funcParams for sorption MT governing equations


                    % add parameters in Freundlich Adsorption model

            end

            % adds a funcParam in the liquid phase and in the sorped phase
            comp.funcParams{end+1} = funcObj;
            comp.sorpFuncParams{end+1} = funcObj;
        end

        % comp: comps class ref, solventName: string
        function comp = removeSorpFuncParamObj(comp, solventName)
            % removes funcParam from both liquid and sorped phase
            for k=1:1:length(comp.sorpFuncParams)
                if strcmp(comp.sorpFuncParams{k}.solventName,solventName)
                    comp.sorpFuncParams(k) = [];
                    break;
                end
            end
            for k=1:1:length(comp.funcParams)
                if isfield(comp.funcParams{k},'solventName')
                    if strcmp(comp.funcParams{k}.solventName,solventName)
                        comp.funcParams(k) = [];
                        break;
                    end
                end
            end
            rmIdx = [];
            for k=1:1:length(comp.sorpHelpers)
                if contains(comp.sorpHelpers{k}.getSubFuncName(),solventName)
                    rmIdx(end+1) = k;
                end
            end
            comp.sorpHelpers(rmIdx) = [];
        end

        % comp: Component class ref, phaseA: string, phaseB: string
        function currGovFuncRaw = getMTGovFunc(comp,phaseA,phaseB)
            if strcmp(phaseA,'Liquid')
                for k=2:1:length(comp.funcParams)
                    if strcmp(comp.funcParams{k}.funcName,[phaseB,' MT'])
                        currGovFuncRaw = comp.funcParams{k}.funcVal;
                    end
                end
            elseif strcmp(phaseA,'Gas')
                for k=1:1:length(comp.gasBulkFuncParams)
                    if strcmp(comp.gasBulkFuncParams{k}.funcName,'Liquid MT')
                        currGovFuncRaw = comp.gasBulkFuncParams{k}.funcVal;
                    end
                end
            else
                for k=1:1:length(comp.sorpFuncParams)
                    if strcmp(comp.sorpFuncParams{k}.solventName,phaseA)
                        currGovFuncRaw = comp.sorpFuncParams{k}.funcVal;
                    end
                end
            end
        end
    end

    methods (Static)
        % num: number, name: string, sym: string,
        % val: number | string, unit: string, funcObj: struct, editable:
        % boolean, defaultParamVals: struct
        function funcObj = createNewParam(locNum, name, sym, val, unit, funcObj, editable, defaultParamVals)
            % won't maintain global number until very end, when all
            % comps and chemicals are set
            paramNames = fields(defaultParamVals);
            if strcmp(sym,"R")
                val = 8.314;
                unit = 'kPa*L/mol*K';
                editable = false;
            elseif strcmp(sym,"tau")
                val = 1;
                unit = 's';
                editable = false;
            elseif strcmp(sym,"pi")
                val = pi;
                unit = '';
                editable = false;
            elseif strcmp(sym,"K_w")
                val = 1E-14;
                unit = 'mol/L';
                editable = false;
            elseif any(strcmp(sym,paramNames))
                val = defaultParamVals.(sym);
                unit = 'g/mol';
                editable = false;
            end

            param = struct( ...
                'locNum',locNum, ...
                'gloNum',1, ...
                'name',name, ...
                'sym',sym, ...
                'val',val, ...
                'unit',unit, ...
                'editable',editable ...
            );
            funcObj.params{end+1} = param;
        end

        % paramSym: string, newVal: number, newUnit: num, newParamName: string, funcObj: struct
        function funcObj = reviseParam(paramSym, newVal, newUnit, newParamName, funcObj)
            paramSyms = cell(1,numel(funcObj.params));
            for k=1:1:numel(funcObj.params)
                paramSyms{k} = funcObj.params{k}.sym;
            end
            idx = strcmp(paramSyms,paramSym);
            funcObj.params{idx}.name = newParamName;
            funcObj.params{idx}.val = newVal;
            funcObj.params{idx}.unit = newUnit;
        end

        % paramName: string, funcObj: struct
        function newFuncObj = removeParams(funcObj)
            funcObj.params = {};
            newFuncObj = funcObj;
        end

        % funcVal: string, funcName: string
        function funcObj = createFuncParamObj(funcVal, funcName)
            funcObj = struct( ...
                'funcVal',funcVal, ...
                'funcName',funcName, ...
                'funcType','None', ...
                'params',struct([]), ...
                'lims',cell(1,1) ...
            );
        end

        % funcVal: string, funcName: string, funcObj: struct
        function newFuncObj = reviseFuncParamObj(funcVal, funcName, funcObj, varargin)
            funcType = '';
            if ~isempty(varargin), funcType = varargin{1}; end
            funcObj.funcVal = funcVal;
            funcObj.funcName = funcName;
            funcObj.funcType = funcType;
            newFuncObj = funcObj;
        end

        % funcParam: funcParamObj, funcVal: string, funcName: string
        function newFuncParam = removeFuncParamObj(funcParam, funcVal, funcName)
            for k=1:1:length(funcParam)
                if strcmp(funcParam{k}.funcVal,funcVal) && ...
                        strcmp(funcParam{k}.funcName,funcName)
                    funcParam(k) = [];
                    newFuncParam = funcParam;
                    return;
                end
            end
        end

        % funcParam: cell{}
        function funcVal = convertPWToLaTeX(funcParam)
            funcVal = funcParam.funcVal;
            lim = funcParam.lims;
            for k=1:1:length(lim)
                newFunc = "\begin{cases} " + lim{k}.funcVal + " & if & " + ...
                    lim{k}.lowerLim + lim{k}.lowLimOp + lim{k}.sysVar + lim{k}.upLimOp + lim{k}.upperLim + ...
                    " \\ " + lim{k}.elseVal + " & else \end{cases}";
                funcVal = replace(funcVal,funcParam.lims{k}.funcVal,newFunc);
            end
        end

        % funcParam: cell{}
        function funcVal = compilePWFuncs(funcParam)
            funcVal = funcParam.funcVal;
            lim = funcParam.lims;
            for k=1:1:length(lim)
                newFunc = "((" + lim{k}.funcVal + ")*(" + lim{k}.lowerLim + ...
                    + lim{k}.lowLimOp + lim{k}.sysVar + "&&" + lim{k}.sysVar + lim{k}.upLimOp + lim{k}.upperLim + ...
                    ")+(" + lim{k}.elseVal + ")*(~(" + lim{k}.lowerLim + ...
                    + lim{k}.lowLimOp + lim{k}.sysVar + "&&" + lim{k}.sysVar + lim{k}.upLimOp + lim{k}.upperLim + ")))";
                funcVal = replace(funcVal,funcParam.lims{k}.funcVal,newFunc);
            end
        end
    end
end