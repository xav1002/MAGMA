classdef Plot < handle
    % Notes:
    % Additional Specifications:
    %   IV: if IC, then need the loEvalLim, upEvalLim, and NbEvalPts
    %   DV: if IC, then need evaltVal, and t is not IV
    %   Units: all variables need unit selection

    properties
        title = ""; % title of plot
        axes = {}; % axes objects

        axesNb = 2; % number of axes on plot

        downloadDir = "None Selected"; % directory to which to download the plot
        fileType = ".png"; % file type the plot should be downloaded as

        display = true; % boolean to determine whether or not to display the plot
        download = false; % boolean to determien whether or not to download the plot
    end

    methods
        % title: string, varNames: string[]
        function plot = Plot(title,varNames)
            plot.title = title;
            for k=1:1:2
                % ### FIXME
                % need to create plot with all standard values in title
                % slots
                if k == 1
                    title = "X Axis Title";
                    isDV = false;
                else
                    title = "Y Axis Title";
                    isDV = true;
                end
                varNames(:,2)
                plot.axes{k} = plot.createNewAxes(k,title,string(varNames(k,2)),string(varNames(:,2)),"min("+varNames(k,2)+")","max("+varNames(k,2)+")", ...
                    true,{0},0,0,50,isDV,false); % ### FIXME
            end
        end

        % plot: Plot class ref, prop: string, val: any
        function updatePlot(plot,prop,val)
            plot.(prop) = val;
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
            for k=1:1:plot.axesNb
                for l=1:1:length(props)
                    axes{k}{l} = plot.getAxProp(k,props(l)); %#ok<AGROW> 
                end
            end
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
            plot.axes{3} = [];
            plot.updateAx(2,"isDV",true);
        end

        % plot: Plot class ref
        function addZAxis(plot)
            varNames = plot.axes{1}.varNameOpts;
            plot.axes{3} = plot.createNewAxes(3,"Z-Axis Title",varNames{3},{varNames},"min("+varNames{3}+")","max("+varNames{3}+")", ...
                true,1,0,0,50,true,false);
            plot.updateAx(2,"isDV",false);
        end

        % plot: Plot class ref, axisName: string
        function axisData = getAxisICEvalData(plot,axisName)
            for k=1:1:length(plot.axes)
                if plot.axes{k}.title == axisName
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
        % loEvalLim: string | number, upEvalLim: string | number, nbEvalPts: number
        function setAxisICEvalData(plot,axisName,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            for k=1:1:length(plot.axes)
                if plot.axes{k}.title == axisName
                    axis = plot.axes{k};
                end
            end
            props = ["evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV"];
            data = {evaltVal,loEvalLim,upEvalLim,nbEvalPts};
            for k=1:1:length(props)
                axis.(props(k)) = data{k};
            end
        end
        
        % plot: Plot class ref, axisName: string, varName: string
        function varData = getVarICEvalData(plot,axisName,varName)
            needDefineEvalt = true;
            for k=1:1:length(plot.axes)
                if plot.axes{k}.title == axisName
                    axis = plot.axes{k};
                end
                if plot.axes{k}.isDV == false && plot.axes{k}.varIsIC == false
                    needDefineEvalt = false;
                end
            end
            for k=1:1:length(axis.varNames)
                if axis.varNames{k} == varName
                    varIdx = k;
                end
            end
            % check if other axes have IC or t in IV
            props = ["evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV","varIsIC"];
            varData = cell(length(props)+1,1);
            for k=1:1:length(varData)
                if k > length(props)
                    varData{k} = needDefineEvalt;
                else
                    if class(axis.(props(k))) == "cell"
                        varData{k} = axis.(props(k)){varIdx};
                    else
                        varData{k} = axis.(props(k));
                    end
                end
            end
        end

        % plot: Plot class ref, axisName: string, varName: string, evaltVal: string | number, 
        % loEvalLim: string | number, upEvalLim: string | number, nbEvalPts: number
        function setVarICEvalData(plot,axisName,varName,evaltVal,loEvalLim,upEvalLim,nbEvalPts)
            for k=1:1:length(plot.axes)
                if plot.axes{k}.title == axisName
                    axis = plot.axes{k};
                end
            end
            for k=1:1:length(axis.varNames)
                if axis.varNames{k} == varName
                    varIdx = k;
                end
            end
            props = ["evaltVal","loEvalLim","upEvalLim","nbEvalPts","isDV","varIsIC"];
            data = {evaltVal,loEvalLim,upEvalLim,nbEvalPts};
            for k=1:1:length(props)
                if class(axis.(props(k))) == "cell"
                    axis.(props(k)){varIdx} = data{k};
                else
                    axis.(props(k)) = data{k};
                end
            end
        end

        % plot: Plot class ref, axes: axes object
        function lims = getDispLims(plot,axes)
            % ### FIXME: code to interpret the lim value from the lim value
            % string
        end
    end

    methods (Static)
        % dir: number, title: string, varName:
        % string, loDispLim: number, upDispLim: number, useDefR: boolean,
        % lowEvalLim: number, upEvalLim: number, nbEvalPts: number, isDV:
        % boolean
        function newAxes = createNewAxes(dir,title,varName,varNameOpts,loDispLim,upDispLim,useDefR,evaltVal,loEvalLim,upEvalLim,nbEvalPts,isDV,varIsIC)
            newAxes = struct( ...
                'dir',dir, ...
                'title',title, ...
                'varNames',varName, ...
                'varNameOpts',varNameOpts, ...
                'loDispLim',loDispLim, ...
                'upDispLim',upDispLim, ...
                'useDefR',useDefR, ...
                'evaltVal',evaltVal, ...
                'loEvalLim',loEvalLim, ...
                'upEvalLim',upEvalLim, ...
                'nbEvalPts',nbEvalPts, ...
                'isDV',isDV, ...
                'varIsIC',varIsIC ...
            );
        end

        function props = getPlotPropNames()
            props = [
                "title", ...
                "axes", ...
                "axesNb", ...
                "downloadDir", ...
                "fileType", ...
                "display", ...
                "download" ...
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
                "evaltVal", ...
                "loEvalLim", ...
                "upEvalLim", ...
                "nbEvalPts", ...
                "isDV", ...
                "varIsIC", ...
            ];
        end
    end
end