classdef Plot < handle
    % Notes:
    % Additional Specifications:
    %   IV: if IC, then need the loEvalLim, upEvalLim, and NbEvalPts
    %   DV: if IC, then need evaltVal, and t is not IV
    %   Units: all variables need unit selection

    % ### FIXME: update default axis titles based on number of dimensions

    properties
        title = ""; % title of plot
        axes = {}; % axes objects
        subplotGroup = 1;
        subplotSlot = 1;

        dimNb = 2; % number of axes on plot

        downloadDir = "None Selected"; % directory to which to download the plot
        fileType = ".png"; % file type the plot should be downloaded as

        display = false; % boolean to determine whether or not to display the plot
        download = false; % boolean to determien whether or not to download the plot

        hold = 'on';
    end

    methods
        % title: string, varNames: string[]
        function plot = Plot(title,varNames,DVNames,group,slot)
            plot.title = title;
            for k=1:1:3
                ICNames = string.empty(0,1);
                for l=1:1:size(DVNames,1)
                    ICNames(l,1) = string(DVNames(l,1))+"~(IC)";
                end
                min = 0;
                max = 100;
                if k == 1
                    title = "X Axis Title";
                    isDV = false;
                    plot.axes{k} = plot.createNewAxes(k,title,string(varNames(k,1)),string(varNames(:,1)),min,max, ...
                        true,{0},{0},{0},{50},isDV,false,ICNames);
                elseif k == 2
                    title = "Y Axis (left) Title";
                    isDV = true;
                    plot.axes{k} = plot.createNewAxes(k,title,string(varNames(k,1)),string(varNames(:,1)),min,max, ...
                        true,{0},{0},{0},{50},isDV,false,ICNames);
                elseif k == 3
                    title = "Y Axis (right) Title";
                    isDV = true;
                    plot.axes{k} = plot.createNewAxes(k,title,"",string(varNames(:,1)),min,max, ...
                        true,{0},{0},{0},{50},isDV,false,ICNames);
                end
            end

            plot.subplotGroup = group;
            plot.subplotSlot = slot;
        end

        % plot: Plot class ref, prop: string, val: any
        function updatePlot(plot,prop,val)
            plot.(prop) = val;
            if strcmp(prop,"dimNb")
                if val == 3
                    plot.axes{2}.varNames = plot.axes{2}.varNames(1);
                    plot.axes{2}.varUnits = plot.axes{2}.varUnits(1);
                end
            end
        end

        % plot: Plot class ref, prop: string
        function val = getPlotProp(plot,prop)
            val = plot.(prop);
        end

        function plotProps = getAllPlotProps(plot)
            plotProps = {};
            props = plot.getPlotPropNames();
            for k=1:1:length(props)
                plotProps{k} = plot.getPlotProp(props(k)); %#ok<AGROW> 
            end
        end

        % plot: Plot class ref, dir: string, prop: string, val: any
        function updateAx(plot,dir,prop,val)
            plot.axes{dir}.(prop) = val;
        end

        % plot: Plot class ref, dir: string, prop: string
        function val = getAxProp(plot,dir,prop)
            val = plot.axes{dir}.(prop);
        end

        % plot: Plot class ref
        function axes = getAllAxProps(plot)
            axes = {};
            props = plot.getAxPropNames();
            for k=1:1:length(plot.axes)
                for l=1:1:length(props)
                    axes{k}{l} = plot.getAxProp(k,props(l)); %#ok<AGROW> 
                end
            end
        end

        % plot: Plot class ref
        function [group,slot] = getSubplotGroupAndSlot(plot)
            group = plot.subplotGroup;
            slot = plot.subplotSlot;
        end

        % plot: Plot class ref
        function downloadPlot(plot,fig)
            % download plot to plot.directory
            try
                saveas(fig,plot.downloadDir+"\"+plot.title+plot.fileType);
            catch err
                throw(err);
            end
        end

        % plot: Plot class ref
        function removeZAxis(plot)
            plot.axes(3) = [];
            plot.updateAx(2,"isDV",true);
        end

        % plot: Plot class ref
        function addZAxis(plot)
            % ### FIXME: add mins and maxes
            varNames = plot.axes{1}.varNameOpts;
            ICNames = plot.axes{1}.ICNames;

            min = 0;
            max = 100;
            plot.axes{3} = plot.createNewAxes(3,"Z-Axis Title",varNames(3),varNames,min,max, ...
                true,1,0,0,50,true,false,ICNames);
            plot.updateAx(2,"isDV",false);
        end

        % % plot: Plot class ref, axisName: string
        function axisData = getAxisICEvalData(plot,axisName)
            for k=1:1:length(plot.axes)
                if strcmp(string(plot.axes{k}.title),string(axisName))
                    axis = plot.axes{k};
                end
            end
            props = ["varNames","evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV","varIsIC"];
            axisData = cell(length(props),1);
            for k=1:1:length(axisData)
                axisData{k} = axis.(props(k));
            end
        end

        % plot: Plot class ref, axisName: string, evaltVal: string | number, 
        % loEvalLim: string | number, upEvalLim: string | number,
        % nbEvalPts: number
        function setAxisICEvalData(plot,axisName,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            % ### FIXME: need to be able to set isDV
            for k=1:1:length(plot.axes)
                if strcmp(string(plot.axes{k}.title),string(axisName))
                    axis = plot.axes{k};
                end
            end
            props = ["evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV","varIsIC"];
            data = {evaltVal,loEvalLim,upEvalLim,nbEvalPts,axis.isDV,axis.varIsIC};
            for k=1:1:length(props)
                axis.(props(k)) = data{k};
            end
        end
        
        % plot: Plot class ref, axisName: string, varName: string
        function varData = getVarICEvalData(plot,axisName,varName)
            needDefineEvalt = true;
            for k=1:1:length(plot.axes)
                if strcmp(string(plot.axes{k}.title),string(axisName))
                    axis = plot.axes{k};
                end
                if ~plot.axes{k}.isDV && ~plot.axes{k}.varIsIC
                    needDefineEvalt = false;
                end
            end
            for k=1:1:length(axis.varNames)
                if strcmp(string(axis.varNames{k}),string(varName))
                    varIdx = k;
                end
            end
            % check if other axes have IC or t in IV
            props = ["varUnits","evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV","varIsIC"];
            varData = cell(length(props)+1,1);
            for k=1:1:length(varData)
                if k > length(props)
                    varData{k} = needDefineEvalt;
                else
                    if iscell(axis.(props(k)))
                        varData{k} = axis.(props(k)){varIdx};
                    else
                        varData{k} = axis.(props(k));
                    end
                end
            end
        end

        % plot: Plot class ref, axisName: string, varName: string, varUnits: string, evaltVal: string | number, 
        % loEvalLim: string | number, upEvalLim: string | number, nbEvalPts: number
        function setVarICEvalData(plot,axisName,varName,varUnits,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            for k=1:1:length(plot.axes)
                if strcmp(string(plot.axes{k}.title),string(axisName))
                    for l=1:1:length(plot.axes{k}.varNames)
                        if strcmp(string(plot.axes{k}.varNames{l}),string(varName))
                            varIdx = l;
                        end
                    end
                    props = ["varUnits","evaltVal","loEvalLim","upEvalLim","nbEvalPts"];
                    data = {varUnits,evaltVal,loEvalLim,upEvalLim,nbEvalPts};
                    for l=1:1:length(props)
                        if iscell(plot.axes{k}.(props(l)))
                            plot.axes{k}.(props(l)){varIdx} = data{l};
                        else
                            plot.axes{k}.(props(l)) = data{l};
                        end
                    end
                    break;
                end
            end
        end

        % plot: Plot class ref, newVarOpts: string[]
        function updateVarOpts(plot,newVarOpts)
            for k=1:1:length(plot.axes)
                plot.axes{k}.varNameOpts = string(newVarOpts(:,1));
            end
        end
    end

    methods (Static)
        % dir: number, title: string, varName:
        % string, loDispLim: number, upDispLim: number, useDefR: boolean,
        % lowEvalLim: number, upEvalLim: number, nbEvalPts: number, isDV:
        % boolean
        function newAxes = createNewAxes(dir,title,varName,varNameOpts,loDispLim,upDispLim,useDefR,evaltVal,loEvalLim,upEvalLim,nbEvalPts,isDV,varIsIC,ICNames)
            newAxes = struct( ...
                'dir',dir, ...
                'title',title, ...
                'varNames',varName, ...
                'varNameOpts',varNameOpts, ...
                'loDispLim',loDispLim, ...
                'upDispLim',upDispLim, ...
                'useDefR',useDefR, ...
                'varUnits',{{''}}, ...
                'evaltVal',{evaltVal}, ...
                'loEvalLim',{loEvalLim}, ...
                'upEvalLim',{upEvalLim}, ...
                'nbEvalPts',{nbEvalPts}, ...
                'isDV',isDV, ...
                'varIsIC',varIsIC, ...
                'ICNames',ICNames ...
            );
        end

        function props = getPlotPropNames()
            props = [
                "title", ...
                "axes", ...
                "subplotGroup", ...
                "subplotSlot", ...
                "dimNb", ...
                "downloadDir", ...
                "fileType", ...
                "display", ...
                "download", ...
                "hold" ...
            ];
        end

        function props = getAxPropNames()
            props = [                
                "dir", ...
                "title", ...
                "varNames", ...
                "varNameOpts", ...
                "loDispLim", ...
                "upDispLim", ...
                "useDefR", ...
                "varUnits", ...
                "evaltVal", ...
                "loEvalLim", ...
                "upEvalLim", ...
                "nbEvalPts", ...
                "isDV", ...
                "varIsIC", ...
                "ICNames" ...
            ];
        end
    end
end