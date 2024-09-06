classdef CompDefaults
    properties (Constant)
        % defaultFuncVals = struct( ...
        %     'Biological Solute', struct( ...
        %                         'Monod', "C_1/(C_1+K_1_1)", ...
        %                         'Linear', "A*C1", ...
        %                         'Custom', "1" ...
        %                     ), ...
        %     'Chemical Solute', struct( ...
        %                         'Monod', "C1/(C1+K1_1)", ...
        %                         'Linear', "A*C1", ...
        %                         'Custom', "1" ...
        %                     ) ...
        % );
    end

    % Format for field - paramVal:
    % replace '^', ' ' with '_'
    % replace '+' with 'p'
    % replace '-' with 'n'
    % replace 'Half-Saturation' with 'K', g/L
    % replace 'Yield' with 'Y', g/g
    % replace 'Death Constant' with 'k', 1/hr
    % replace 'Inhibition Constant' with 'KI', g/L

    methods (Static)
        function specD = CompDefaults()
%             specD.Synechoccocus_elongatus_UTEX_2973.K_NO3_n = 0.032;
        end

        function [expr,descr,params] = getDefaultFuncVal(comp,funcType)
            if strcmp(comp.getType(),'Biological Solute')
                bioName = comp.getName();
                bioNum = comp.getNum();
                chemName = "Substrate";
                chemNum = 1;
            elseif strcmp(comp.getType(),'Chemical Solute')
                chemName = comp.getName();
                chemNum = comp.getNum();
                bioName = "Biological Component";
                chemNum = 1;
            end
            switch funcType
                case 'Monod'
                    expr = char("mu_max_" + bioNum + "*(C_"+chemNum+"/(K_"+bioNum+"_"+chemNum+"+C_"+chemNum+"))");
                    descr = "Monod model expression that relates the concentration of $\mathit{\mathbf{C_"+chemNum+"}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                case 'Moser'
                    expr = char("mu_max_" + bioNum + "*((C_"+chemNum+"^(n_"+bioNum+"_"+chemNum+"))/(K_"+bioNum+"_"+chemNum+"+C_"+chemNum+"^(n_"+bioNum+"_"+chemNum+")))");
                    descr = "Moser model expression that relates the concentration of $\mathit{\mathbf{C_"+chemNum+"}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_[# chemical substrate] $\\$ n_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{"+chemNum+"}}","n_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(3))) + "$: Exponential parameter.";
                case 'Contois'
                    expr = char("mu_max_" + bioNum + "*(C_"+chemNum+"/(K_" + bioNum + "_"+chemNum+"*X_"+bioNum+"+C_"+chemNum+"))");
                    descr = "Contois model expression that relates the concentration of $\mathit{\mathbf{C_"+chemNum+"}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Contois model constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                case 'Haldane'
                    expr = char("mu_max_" + bioNum + "*(C_"+chemNum+"/(K_"+bioNum+"_"+chemNum+"+C_"+chemNum+"))*(K_i_"+bioNum+"_"+chemNum+"/(K_i_"+bioNum+"_"+chemNum+"+C_"+chemNum+"))");
                    descr = "Haldane model expression that relates the concentration of $\mathit{\mathbf{C_"+chemNum+"}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$ while considering substrate inhibition of $\mathit{\mathbf{C_"+chemNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_[# chemical substrate] $\\$ K_i_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{"+chemNum+"}}","K_i_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(3))) + "$: Inhibition constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                case 'Grant'
                    expr = char("mu_max_" + bioNum + "*(1/(K_i_"+bioNum+"_"+chemNum+"+C_"+chemNum+"))");
                    descr = "Grant model expression that considers the inhibition effect of $\mathit{\mathbf{C_"+chemNum+"}}$ on the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_i_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_i_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Inhibition constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                case 'Andrews'
                    expr = char("mu_max_" + bioNum + "*(C_"+chemNum+"/(K_"+bioNum+"_"+chemNum+"+C_"+chemNum+"+(C_"+chemNum+"^2)/(K_i_"+bioNum+"_"+chemNum+")))");
                    descr = "Andrews model expression that relates the concentration of $\mathit{\mathbf{C_"+chemNum+"}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$ while considering substrate inhibition of $\mathit{\mathbf{C_"+chemNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_[# chemical substrate] $\\$ K_i_[# biological component]_[# chemical substrate]";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{"+chemNum+"}}","K_i_{"+bioNum+"_{"+chemNum+"}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(3))) + "$: Inhibition constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{C_"+chemNum+"}}$";
                case 'pH Model'
                    expr = char("mu_max_" + bioNum + "/(1+K_i_"+bioNum+"_H3O/H3O+K_"+bioNum+"_H3O*H3O)");
                    descr = "pH model expression that relates the concentration of $\mathit{\mathbf{H_3O^+}}$ to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$ while considering inhibition effect of $\mathit{\mathbf{H_3O^+}}$ on growth rate. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_i_[# biological component]_H3O $\\$ K_[# biological component]_H3O";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_{H3O}}","K_i_{"+bioNum+"_{H3O}}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{H_3O^+}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(3))) + "$: Inhibition constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to $\mathit{\mathbf{H_3O^+}}$";
                % case 'pH Model 2'
                % case 'pH Model 3'
                case 'Light Intensity: Tamiya'
                    expr = char("mu_max_" + bioNum + "*(I/(K_"+bioNum+"_I+I))");
                    descr = "Tamiya model expression that relates \textbf{average light intensity} to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_I";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_I}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to \textbf{average light intensity}";
                case 'Light Intensity: van Oorschot'
                    expr = char("mu_max_" + bioNum + "*exp(1-I/K_"+bioNum+"_I)");
                    descr = "van Oorschot model expression that relates \textbf{average light intensity} to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_I";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_I}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: van Oorschot model constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to \textbf{average light intensity}";
                % case 'Light Intensity: Chalker'
                case 'Light Intensity: Aiba'
                    expr = char("mu_max_" + bioNum + "*(I/(K_"+bioNum+"_I+I+(I^2)/(K_i_"+bioNum+"_I)))");
                    descr = "Tamiya model expression that relates \textbf{average light intensity} to the specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$. $\\$" + ...
                        "Recommended Notation: $\\$ $\mu$_max_[# biological component] $\\$ K_[# biological component]_I";

                    param_syms = ["mu_max_{"+bioNum+"}","K_{"+bioNum+"_I}","K_i_{"+bioNum+"_I}"];
                    params = "$" + uni2latex(char(param_syms(1))) + "$: Maximum specific growth rate of $\mathit{\mathbf{X_"+bioNum+"}}$";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Half-saturation constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to \textbf{average light intensity}";
                    params = params + "$\\" + uni2latex(char(param_syms(2))) + "$: Inhibition constant for $\mathit{\mathbf{X_"+bioNum+"}}$ growth rate with respect to \textbf{average light intensity}";
            end
        end

        function MTFuncVal = getDefaultMTFuncVals(MTFuncName,phaseACompSym,phaseBCompSym,phaseASym,phaseBSym,compMWSym)
            switch MTFuncName
                case "None"
                    MTFuncVal = "0";

                case "Inter-Phase Transfer based on Liquid Phase Overall Mass Transfer Coefficient and Driving Force"
                    compNum = regexp(phaseACompSym,'\d*','match');
                    phaseLetter = ["S","S"];
                    if strcmp(phaseASym,"V")
                        phaseLetter(1) = "L";
                    end
                    if strcmp(phaseBSym,"V")
                        phaseLetter(2) = "L";
                    end
                    if strcmp(phaseASym,"(V_m-V_tot)")
                        phaseLetter(1) = "G";
                    end
                    if strcmp(phaseBSym,"(V_m-V_tot)")
                        phaseLetter(2) = "G";
                    end

                    % ### FIXME: need to make sure that the mass transfer
                    % coefficient is the same for both phases
                    if strcmp(phaseASym,"(V_m-V_tot)")
                        MTFuncVal = "k_"+phaseLetter(1)+"_"+phaseLetter(2)+"_"+compNum{1}+"*("+"R*T"+")/("+phaseASym+"*"+compMWSym+")*("+phaseBCompSym+"-"+phaseACompSym+")";
                    elseif strcmp(phaseASym,"V")
                        MTFuncVal = "k_"+phaseLetter(1)+"_"+phaseLetter(2)+"_"+compNum{1}+"/"+phaseASym+"*("+phaseBCompSym+"-"+phaseACompSym+")";
                    else
                        MTFuncVal = "k_"+phaseLetter(1)+"_"+phaseLetter(2)+"_"+compNum{1}+"/"+phaseASym+"*("+phaseBCompSym+"-"+phaseACompSym+")";
                    end

                case "Inter-Phase Transfer based on Liquid Phase Overall Capacity Coefficient and Driving Force"
                    compNum = regexp(phaseACompSym,'\d*','match');
                    phaseLetter = ["S","S"];
                    if strcmp(phaseASym,"V")
                        phaseLetter(1) = "L";
                    end
                    if strcmp(phaseBSym,"V")
                        phaseLetter(2) = "L";
                    end
                    if strcmp(phaseASym,"(V_m-V_tot)")
                        phaseLetter(1) = "G";
                    end
                    if strcmp(phaseBSym,"(V_m-V_tot)")
                        phaseLetter(2) = "G";
                    end

                    if strcmp(phaseASym,"(V_m-V_tot)")
                        MTFuncVal = "0";
                    elseif strcmp(phaseASym,"V")
                        MTFuncVal = "K_"+phaseLetter(1)+"_"+phaseLetter(2)+"_"+compNum{1}+"*("+phaseBCompSym+"-"+phaseACompSym+")";
                    else
                        MTFuncVal = "K_"+phaseLetter(1)+"_"+phaseLetter(2)+"_"+compNum{1}+"/"+phaseASym+"*("+phaseBCompSym+"-"+phaseACompSym+")";
                    end

            end
        end
    end
end