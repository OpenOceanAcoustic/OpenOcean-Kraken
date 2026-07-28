function tests = test_m_kraken
tests = functiontests(localfunctions);
end

function testBiologicalAttenuationJsonAndCli(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'MunkK.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>

baseRoot = fullfile(work, 'base');
command = sprintf('"%s" "%s" --output "%s" --threads 1 --mod-only --json', ...
    exe, fixture, baseRoot);
[status, output] = system(command);
verifyEqual(testCase, status, 0, output);

model = krakenDataModel(exe);
model.loadJson([baseRoot '.json']);
model.setBiologicalAttenuation([
    10.0 30.0 1000.0 5.0 0.04
    40.0 60.0 1200.0 4.0 0.02
]);
jsonPath = model.Write(fullfile(work, 'biological.json'));
payload = jsondecode(fileread(jsonPath));
verifyEqual(testCase, payload.AttenUnit.AttenuationUnit, 'dB/lambda');
verifyEqual(testCase, payload.AttenUnit.OceanAbsorptionModel, 'Biological');
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).Z1, 10.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).Z2, 30.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).f0, 1000.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).Q, 5.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(1).a0, 0.04);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).Z1, 40.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).Z2, 60.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).f0, 1200.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).Q, 4.0);
verifyEqual(testCase, payload.AttenUnit.BiologicalLayers(2).a0, 0.02);

resultRoot = fullfile(work, 'validated');
command = sprintf('"%s" "%s" --output "%s" --threads 1 --mod-only', ...
    exe, jsonPath, resultRoot);
[status, output] = system(command);
verifyEqual(testCase, status, 0, output);
verifyTrue(testCase, isfile([resultRoot '.mod']));

model.setBiologicalAttenuation(zeros(0, 5));
emptyJson = model.Write(fullfile(work, 'biological_empty.json'));
emptyPayload = jsondecode(fileread(emptyJson));
verifyTrue(testCase, isempty(emptyPayload.AttenUnit.BiologicalLayers));

verifyError(testCase, ...
    @() model.setBiologicalAttenuation([30 10 1000 5 0.04]), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation([10 30 1000 0 0.04]), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation([]), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation({}), ...
    'OpenOceanKraken:InvalidBiologicalLayer');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), 'dB/l'), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), 'db/lambda'), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), []), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), ''), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), string.empty), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
verifyError(testCase, ...
    @() model.setBiologicalAttenuation(zeros(0, 5), {}), ...
    'OpenOceanKraken:InvalidAttenuationUnit');
end

function testEnvCliAndReaders(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'MunkK.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>

model = krakenDataModel(exe);
model.NumThreads = 1;
model.loadEnv(fixture);
resultRoot = model.run(work, 'Velocity', true, 'Mod', true);
pressure = model.getPressure(resultRoot);
vertical = model.getVerticalVelocity(resultRoot);
modes = model.getModes(resultRoot);

verifyEqual(testCase, ndims(pressure.values), 3);
verifyEqual(testCase, size(vertical.values), size(pressure.values));
verifyTrue(testCase, all(isfinite(real(pressure.values)), 'all'));
verifyGreaterThan(testCase, numel(modes.profiles), 0);
verifyGreaterThan(testCase, numel(modes.profiles(1).wavenumbers), 0);
end

function testJsonRoundTripAndPlots(testCase)
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
exe = fullfile(root, 'bin', 'OpenOceanKraken.exe');
fixture = fullfile(fileparts(root), 'test', 'MunkK.env');
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>

model = krakenDataModel(exe);
model.NumThreads = 1;
model.loadEnv(fixture);
resultRoot = model.run(work, 'Json', true, 'Mod', true);

copy = krakenDataModel(exe);
copy.NumThreads = 1;
copy.loadJson([resultRoot '.json']);
jsonPath = copy.Write(work);
verifyTrue(testCase, isfile(jsonPath));

pressureFigure = figure('Visible', 'off');
modeFigure = figure('Visible', 'off');
figureCleanup = onCleanup(@() close([pressureFigure, modeFigure])); %#ok<NASGU>
model.plotPressure(resultRoot, 'Figure', pressureFigure);
model.plotModes(resultRoot, 'Figure', modeFigure);
verifyGreaterThan(testCase, numel(findall(pressureFigure, 'Type', 'axes')), 0);
verifyGreaterThan(testCase, numel(findall(modeFigure, 'Type', 'axes')), 1);
end

function testTruncatedReaders(testCase)
work = tempname;
mkdir(work);
cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
shd = fullfile(work, 'bad.shd');
mod = fullfile(work, 'bad.mod');
writeBytes(shd, zeros(1, 32, 'uint8'));
writeBytes(mod, zeros(1, 32, 'uint8'));
verifyError(testCase, @() read_shd(shd), 'OpenOceanKraken:InvalidSHD');
verifyError(testCase, @() read_mod(mod), 'OpenOceanKraken:InvalidMOD');
end

function writeBytes(path, values)
fid = fopen(path, 'wb');
assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, values, 'uint8');
end
