% h = erp_plot_signal_fromclass(class, epochs, varargin)
%
%       Plots an ERP class activation signal. In your colour map first's 
%       colour, solid line: the mean signal as defined. The dotted and
%       dashed lines indicate the signal's variability as per the defined
%       deviations. If a slope has been defined, the second colour curves
%       indicate the signal (mean and extremes) at the end of the slope
%       (i.e. the final epoch).
%
% In:
%       class - 1x1 struct, the class variable
%       epochs - single epoch configuration struct containing at least
%                sampling rate in Hz (field name 'srate'), epoch length in ms
%                 ('length'), and the total number of epochs ('n')
%
% Optional (key-value pairs):
%       newfig - (0|1) whether or not to open a new figure window.
%                default: 1
%       baseonly - (0|1) whether or not to plot only the base signal,
%                  without any deviations or sloping. default: 0
%       shadeplot - (0|1) if ploting also deviations and sloping, plot as a
%                   shaded area instead of other colored lines (default:0)
%
% Out:  
%       h - handle of the generated figure
%
% Usage example:
%       >> epochs = struct('n', 100, 'srate', 1000, 'length', 1000);
%       >> erp = struct('type', 'erp', 'peakLatency', 300, ...
%       >>      'peakLatencySlope', 150, 'peakWidth', 300, ...
%       >>      'peakWidthDv', 50, 'peakWidthSlope', 100, ...
%       >>      'peakAmplitude', 1, 'peakAmplitudeDv', .25, ...
%       >>      'peakamplitudeSlope', -.25);
%       >> plot_signal_fromclass(erp, epochs);
% 
%                    Copyright 2017, 2018 Laurens R Krol
%                    Team PhyPA, Biological Psychology and Neuroergonomics,
%                    Berlin Institute of Technology

% 2018-06-04 lrk
%   - Added peakLatencyShift
% 2017-06-13 First version

% This file is part of Simulating Event-Related EEG Activity (SEREEGA).

% SEREEGA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% SEREEGA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with SEREEGA.  If not, see <http://www.gnu.org/licenses/>.

function h = erp_plot_signal_fromclass(class, epochs, varargin)

% parsing input
p = inputParser;

addRequired(p, 'class', @isstruct);
addRequired(p, 'epochs', @isstruct);

addParameter(p, 'newfig', 1, @isnumeric);
addParameter(p, 'baseonly', 0, @isnumeric);
addParameter(p, 'shadeplot', 0, @isnumeric);

parse(p, class, epochs, varargin{:})

class = p.Results.class;
epochs = p.Results.epochs;
newfig = p.Results.newfig;
baseonly = p.Results.baseonly;
shadeplot = p.Results.shadeplot;

% getting time stamps
x = 1:epochs.length/1000*epochs.srate;
x = x/epochs.srate;

% correcting time stamps if a prestimulus period is indicated
if isfield(epochs, 'prestim')
    x = x - epochs.prestim/1000; end

if newfig, h = figure('name', 'ERP signal', 'NumberTitle', 'off', 'ToolBar', 'none'); else, h = gcf; end
hold on;

if ~baseonly
    if shadeplot
        m = NaN(1,epochs.length/1000*epochs.srate);M = NaN(1,epochs.length/1000*epochs.srate);
        nsim = 100;
        for i = 1:nsim * epochs.n
            tmp = erp_generate_signal_fromclass(class,epochs,'baseonly',0,'epochNumber',ceil(i/nsim));
            m = min(m,tmp);
            M = max(M,tmp);
        end
        shadebetween(x,M,m,[.8 .8 .8],'none',.5);
    else
        % signal with maximum possible deviation (negative)
        ax = gca;
        colidx = ax.ColorOrderIndex;
        ax.ColorOrderIndex = colidx;
        plot(x, erp_generate_signal( ...
            class.peakLatency - class.peakLatencyDv - class.peakLatencyShift, ...
            class.peakWidth - class.peakWidthDv, ...
            class.peakAmplitude - class.peakAmplitudeDv, ...
            epochs.srate, epochs.length), ':');
        
        % signal with maximum possible deviation (positive)
        ax.ColorOrderIndex = colidx;
        plot(x, erp_generate_signal( ...
            class.peakLatency + class.peakLatencyDv + class.peakLatencyShift, ...
            class.peakWidth + class.peakWidthDv, ...
            class.peakAmplitude + class.peakAmplitudeDv, ...
            epochs.srate, epochs.length), '--');
        
        if any([class.peakLatencySlope, class.peakWidthSlope, class.peakAmplitudeSlope])
            % additionally plotting the signal with maximum slope applied
            % mean
            plot(x, erp_generate_signal(...
                class.peakLatency + class.peakLatencySlope, ...
                class.peakWidth + class.peakWidthSlope, ...
                class.peakAmplitude + class.peakAmplitudeSlope, ...
                epochs.srate, epochs.length), '-');
            
            % negative deviation
            ax.ColorOrderIndex = 2;
            plot(x, erp_generate_signal( ...
                class.peakLatency + class.peakLatencySlope - class.peakLatencyDv, ...
                class.peakWidth + class.peakWidthSlope - class.peakWidthDv, ...
                class.peakAmplitude + class.peakAmplitudeSlope - class.peakAmplitudeDv, ...
                epochs.srate, epochs.length), ':');
            
            % positive deviation
            ax.ColorOrderIndex = 2;
            plot(x, erp_generate_signal( ...
                class.peakLatency + class.peakLatencySlope + class.peakLatencyDv, ...
                class.peakWidth + class.peakWidthSlope + class.peakWidthDv, ...
                class.peakAmplitude + class.peakAmplitudeSlope + class.peakAmplitudeDv, ...
                epochs.srate, epochs.length), '--');
        end
    end
end

% plotting the mean signal, no deviations applied
plot(x, erp_generate_signal_fromclass(class, epochs, 'baseonly', 1), '-');

end

function[fillhandle] = shadebetween(xpoints,upper,lower,color,edge,transparency)
%USAGE: [fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,transparency)
%This function will fill a region with a color between the two vectors provided
%using the Matlab fill command.
%
%fillhandle is the returned handle to the filled region in the plot.
%xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
%upper = the upper curve values (data can be less than lower)
%lower = the lower curve values (data can be more than upper)
%color = the color of the filled area 
%edge  = the color around the edge of the filled area
%transparency = value ranging from 1 for opaque to 0 for invisible for
%       the filled color only.
%
%John A. Bockstege November 2006;
%Example:
%     a=rand(1,20);%Vector of random data
%     b=a+2*rand(1,20);%2nd vector of data points;
%     x=1:20;%horizontal vector
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
%     grid on
%     legend('Datr')
if not(exist('transparency','var')) || isempty(transparency)
    transparency=.5;
end
if not(exist('edge','var')) || isempty(edge)
    edge = 'k';
end
if not(exist('color','var')) || isempty(color)
    color = 'b';
end
if isvector(xpoints)
    xpoints = xpoints(:);
end
if isvector(upper)
    upper = upper(:);
end
if isvector(lower)
    lower = lower(:);
end
filled=[upper;flipud(lower)];
xpoints=[xpoints;flipud(xpoints)];
nans = isnan(filled) | isnan(xpoints);
filled(nans) = [];
xpoints(nans) = [];

fillhandle=fill(xpoints,filled,color,...
    'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);

if nargout == 0
    clear fillhandle
end
end