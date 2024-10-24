classdef Stream < handle
    properties
        name = '';
        phase = '';
        dir = '';
        num = 0;

        sym = '';
        
        compsData = {};
        basis = '';
        basisU = '';
        solventName = '';
        flowrate = 0;
        flowrateU = '';
    end

    methods
        function strm = Stream(name,phase,dir,num,comps)
            if isempty(name), strm.name = ['Stream ',char(string(num))]; else, strm.name = name; end
            strm.phase = phase;
            strm.dir = dir;
            strm.num = num;
            strm.sym = ['STRM_',char(string(num))];

            if strcmp(phase,'L')
                strm.basis = 'Concentration';
                strm.basisU = 'g/L';
                strm.flowrateU = 'L/s';
                strm.solventName = 'Water';
            elseif strcmp(phase,'G')
                strm.basis = 'Mass Fraction';
                strm.basisU = '';
                strm.solventName = '';
                if strcmp(dir,'in')
                    strm.flowrateU = 'g/s';
                elseif strcmp(dir,'out')
                    strm.flowrateU = 'bar/s';
                end
            end

            if strcmp(dir,'in')
                strm.updateAvailComps(comps);
            end
        end

        function name = getName(strm)
            name = strm.name;
        end

        function phase = getPhase(strm)
            phase = strm.phase;
        end

        function dir = getDir(strm)
            dir = strm.dir;
        end

        function num = getNum(strm)
            num = strm.num;
        end

        function sym = getSym(strm)
            sym = strm.sym;
        end

        function [compsData,basis,basisU,solventName,flowrate,flowrateU] = getCompsData(strm)
            compsData = strm.compsData;
            basis = strm.basis;
            basisU = strm.basisU;
            solventName = strm.solventName;
            flowrate = strm.flowrate;
            flowrateU = strm.flowrateU;
        end

        function setName(strm,name)
            strm.name = name;
        end

        function setPhase(strm,phase)
            strm.phase = phase;
        end

        function setDir(strm,dir)
            strm.dir = dir;
        end

        function setNum(strm,num)
            strm.num = num;
        end

        function setSym(strm,sym)
            strm.sym = sym;
        end

        function setCompsData(strm,compsData,basis,basisU,solventName,flowrate,flowrateU)
            strm.compsData = compsData;
            strm.basis = basis;
            strm.basisU = basisU;
            strm.solventName = solventName;
            strm.flowrate = flowrate;
            strm.flowrateU = flowrateU;
        end

        function updateAvailComps(strm,comps)
            if strcmp(strm.phase,'G')
                newComps = {};
                for k=1:1:length(comps)
                    if comps{k}.is_vol
                        newComps{end+1} = comps{k}; %#ok<AGROW>
                    end
                end
                comps = newComps;
            end
            newData = cell(length(comps),2);
            for k=1:1:length(comps)
                newData{k,1} = comps{k}.getName();
                if ~isempty(strm.compsData)
                    prev_val_idx = strcmp(strm.compsData(:,1),newData{k,1});
                else
                    prev_val_idx = 0;
                end
                if any(prev_val_idx)
                    newData{k,2} = strm.compsData{prev_val_idx,2};
                else
                    newData{k,2} = 0;
                end
            end
            strm.compsData = newData;
        end
    end
end