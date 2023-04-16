classdef Chemical < handle
    properties
        number = 0; % chem identifying number
        name = ""; % chem name

        transOut = {}; % struct representing mass transport of Chemical out of system
        initConc = 1; % initial concentration for Chemical

        funcParams = {}; % cell array of structs that contain the values, names, symbols, and units of parameters,
                        % with the indicies being the parameter numbers.

        sym = "";
        

        %% Chemical Struct Data Structure
        % chem:
        %   number: number
        %   name: name
        %   modelType: number that corresponds to the specified model
        %       function
        %   coeffs: array of coefficients required for the specified model
        %       function

        %% Reaction Struct Data Structure
        % rxn:
        %   params: number[]
        %   modelType: number that corresponds to the specified reaction
        %       model

        %% Transport Struct Data Structure
        % trans:
        %   params: number[]
        %   modelType: number that corresponds to the specified transport
        %       model

        rateEqn = 1; % overall dC/dt function for Chemical, all relevant transport and reaction terms will be compiled into this function

        %% Chemical Name Expression Replacements
        % replace ' ' with '_'
        % replace '^', ' ' with '_'
        % replace '+' with 'p'
        % replace '-' with 'n'
    end

    methods
        function chem = Chemical(name, number, initConc)
            chem.name = name;
            chem.number = number;
            chem.initConc = initConc;
            chem.sym = "X"+number;
            
            chem.setDefaultParams(name);
        end

        % chem: Chemical class ref, name: string
        function setDefaultParams(chem,name)

        end

        % chem: Chemical class ref, conc: number
        function chem = setInitConc(chem, conc)
            chem.initConc = conc;
        end

        % chem: Chemical class ref, chemName: string, interact: string,
        % func: function_handle, varNames: string[], paramNames: string[],
        % funcNames: string, funcType: string
        function chem = setModel(chem, chemName, interact, func, varNames, paramNames, funcName, funcType)
            
        end

        function funcNames = getFuncNames(chem)
            funcNames = cell(1,length(chem.funcParams));
            for k=1:1:length(chem.funcParams)
                funcNames{k} = chem.funcParams{k}.funcName;
            end
        end
    end

end