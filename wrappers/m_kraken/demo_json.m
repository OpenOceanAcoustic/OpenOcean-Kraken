function resultRoot = demo_json(outputDir)
%DEMO_JSON Produce a normalized OOK JSON input, reload it, and run it.

wrapperDir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(wrapperDir));
if nargin < 1
    outputDir = fullfile(root, 'tmp', 'm_kraken_demo_json');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end
executable = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'MunkK.env');

% Let OOK normalize the legacy ENV once; the second run is JSON-driven.
bootstrap = krakenDataModel(executable);
bootstrap.NumThreads = 1;
bootstrap.loadEnv(fixture);
bootstrapRoot = bootstrap.run(fullfile(outputDir, 'bootstrap'), 'Json', true);

model = krakenDataModel(executable);
model.NumThreads = 1;
model.loadJson([bootstrapRoot '.json']);
jsonPath = model.Write(outputDir);
model.loadJson(jsonPath);
resultRoot = model.run(fullfile(outputDir, 'json_run'), 'Mod', true);

pressureFigure = figure('Visible', 'off');
modeFigure = figure('Visible', 'off');
cleanup = onCleanup(@() close([pressureFigure, modeFigure])); %#ok<NASGU>
model.plotPressure(resultRoot, 'Figure', pressureFigure);
model.plotModes(resultRoot, 'Figure', modeFigure);
exportgraphics(pressureFigure, fullfile(outputDir, 'pressure_json.png'));
exportgraphics(modeFigure, fullfile(outputDir, 'modes_json.png'));
fprintf('OOK JSON results written under %s\n', resultRoot);
end
