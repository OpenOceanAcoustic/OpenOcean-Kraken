function tests = test_m_kraken
tests = functiontests(localfunctions);
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
