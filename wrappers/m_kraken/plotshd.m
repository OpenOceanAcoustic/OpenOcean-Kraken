function axesHandle = plotshd(field, varargin)
%PLOTSHD Plot an OOK field as transmission loss.

parser = inputParser;
addParameter(parser, 'Figure', [], @(x) isempty(x) || isgraphics(x, 'figure'));
addParameter(parser, 'SourceIndex', 1, @(x) isscalar(x) && x >= 1 && x == floor(x));
addParameter(parser, 'Component', 'pressure', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
sourceIndex = parser.Results.SourceIndex;
if sourceIndex > size(field.values, 1)
    error('OpenOceanKraken:InvalidPlot', 'SourceIndex is outside the field.');
end
if isempty(parser.Results.Figure)
    figureHandle = figure;
else
    figureHandle = parser.Results.Figure;
end
axesHandle = axes('Parent', figureHandle);
value = squeeze(field.values(sourceIndex, :, :));
transmissionLoss = -20 * log10(max(abs(value), realmin('single')));
if startsWith(lower(field.plot_type), 'irregular')
    scatter(axesHandle, field.receiver_ranges / 1000, field.receiver_depths, ...
        10, transmissionLoss(:, 1), 'filled');
else
    surf(axesHandle, field.receiver_ranges / 1000, field.receiver_depths, ...
        transmissionLoss.', 'EdgeColor', 'none');
    view(axesHandle, 2);
end
xlabel(axesHandle, 'Range (km)');
ylabel(axesHandle, 'Depth (m)');
set(axesHandle, 'YDir', 'reverse');
title(axesHandle, sprintf('%s - %.6g Hz', field.title, field.frequency));
component = lower(char(parser.Results.Component));
switch component
    case 'pressure'
        colorLabel = 'Transmission loss (dB)';
    case 'vertical_velocity'
        colorLabel = 'Vertical-velocity level (dB)';
    case 'horizontal_velocity'
        colorLabel = 'Horizontal-velocity level (dB)';
    otherwise
        error('OpenOceanKraken:InvalidPlot', 'Unknown field component.');
end
bar = colorbar(axesHandle);
bar.Label.String = colorLabel;
end
