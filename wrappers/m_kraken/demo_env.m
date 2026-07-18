function resultRoot = demo_env(outputDir)
%DEMO_ENV Run the checked-in ENV example through the MATLAB wrapper.

wrapperDir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(wrapperDir));
if nargin < 1
    outputDir = fullfile(root, 'tmp', 'm_kraken_demo_env');
end
model = krakenDataModel(fullfile(root, 'bin', 'OpenOceanKraken.exe'));
model.NumThreads = 1;
model.loadEnv(fullfile(fileparts(root), 'test', 'MunkK.env'));
resultRoot = model.run(outputDir, 'Mod', true);

pressureFigure = figure('Visible', 'off');
modeFigure = figure('Visible', 'off');
cleanup = onCleanup(@() close([pressureFigure, modeFigure])); %#ok<NASGU>
model.plotPressure(resultRoot, 'Figure', pressureFigure);
model.plotModes(resultRoot, 'Figure', modeFigure);
exportgraphics(pressureFigure, fullfile(outputDir, 'pressure.png'));
exportgraphics(modeFigure, fullfile(outputDir, 'modes.png'));
fprintf('OOK ENV results written under %s\n', resultRoot);
end
