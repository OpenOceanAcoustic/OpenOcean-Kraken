function figureHandle = plotmode(modes, varargin)
%PLOTMODE Plot wavenumbers and selected normalized mode shapes.

parser = inputParser;
addParameter(parser, 'Figure', [], @(x) isempty(x) || isgraphics(x, 'figure'));
addParameter(parser, 'ProfileIndex', 1, @(x) isscalar(x) && x >= 1 && x == floor(x));
addParameter(parser, 'ModeIndices', [], @isnumeric);
parse(parser, varargin{:});
profileIndex = parser.Results.ProfileIndex;
if profileIndex > numel(modes.profiles)
    error('OpenOceanKraken:InvalidPlot', 'ProfileIndex is outside the mode data.');
end
profile = modes.profiles(profileIndex);
if isempty(parser.Results.ModeIndices)
    selected = 1:min(6, numel(profile.wavenumbers));
else
    selected = parser.Results.ModeIndices;
end
if isempty(selected) || any(selected < 1) || any(selected > numel(profile.wavenumbers))
    error('OpenOceanKraken:InvalidPlot', 'ModeIndices contain an invalid mode.');
end
if isempty(parser.Results.Figure)
    figureHandle = figure;
else
    figureHandle = parser.Results.Figure;
end

waveAxes = subplot(1, 2, 1, 'Parent', figureHandle);
modeNumbers = 1:numel(profile.wavenumbers);
plot(waveAxes, modeNumbers, real(profile.wavenumbers), '-', ...
    modeNumbers, imag(profile.wavenumbers), '--');
xlabel(waveAxes, 'Mode number');
ylabel(waveAxes, 'Wavenumber (1/m)');
legend(waveAxes, {'Re(k)', 'Im(k)'}, 'Location', 'best');
grid(waveAxes, 'on');

shapeAxes = subplot(1, 2, 2, 'Parent', figureHandle);
hold(shapeAxes, 'on');
for index = selected
    shape = profile.mode_shapes(index, :);
    scale = max(abs(shape));
    if scale > 0
        shape = real(shape) / scale;
    else
        shape = real(shape);
    end
    plot(shapeAxes, shape, profile.depth, 'DisplayName', sprintf('Mode %d', index));
end
hold(shapeAxes, 'off');
xlabel(shapeAxes, 'Normalized real mode shape');
ylabel(shapeAxes, 'Depth (m)');
set(shapeAxes, 'YDir', 'reverse');
legend(shapeAxes, 'Location', 'best');
grid(shapeAxes, 'on');
sgtitle(figureHandle, sprintf('%s - %.6g Hz', modes.title, modes.frequency));
end
