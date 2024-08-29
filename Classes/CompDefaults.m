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

        function defaultFuncValue = getDefaultFuncVals(comp,funcType)
            compType = comp.getType();
            switch compType
                case 'Biological Solute'
                    switch funcType
                        case 'Monod'
                    end
                case 'Chemical Solute'
                    switch funcType

                    end
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