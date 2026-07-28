classdef krakenDataModel < handle
    %KRAKENDATAMODEL File-driven MATLAB interface for OpenOceanKraken.

    properties
        Title = 'OpenOceanKraken case'
        NumThreads = 1
        Frequency = 0
        FrequencyVector = []
        SourceDepth = []
        ReceiverDepth = []
        ReceiverRange = []
        ProfileRange = 0
        PhaseSpeed = [0, 0]
        SSP = []
        Attenuation = struct()
        GridType = 'Rectangular'
        SourceType = 'Point'
        RunMode = 'Both'
        CoherenceType = 'Coherent'
        ModeType = 'Adiabatic'
        MLimit = 9999
        VelocityEnabled = false
        ReflectionCoef = struct()
        SBP = struct()
        Rmax = 0
    end

    properties (SetAccess = private)
        Executable
    end

    properties (Access = private)
        InputPath = ''
        RawConfig = struct()
        LastResultRoot = ''
    end

    methods
        function obj = krakenDataModel(executable)
            if nargin < 1 || isempty(executable)
                root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
                executable = fullfile(root, 'bin', 'OpenOceanKraken.exe');
            end
            executable = char(executable);
            if ~isfile(executable)
                error('OpenOceanKraken:MissingExecutable', ...
                    'OpenOceanKraken executable does not exist: %s', executable);
            end
            obj.Executable = executable;
        end

        function loadEnv(obj, path)
            path = char(path);
            if ~isfile(path)
                error('OpenOceanKraken:MissingInput', 'ENV input does not exist: %s', path);
            end
            [~, ~, extension] = fileparts(path);
            if ~strcmpi(extension, '.env')
                error('OpenOceanKraken:InvalidInput', 'Expected an ENV input file.');
            end
            obj.InputPath = path;
            obj.RawConfig = struct();
        end

        function loadJson(obj, path)
            path = char(path);
            if ~isfile(path)
                error('OpenOceanKraken:MissingInput', 'JSON input does not exist: %s', path);
            end
            try
                config = jsondecode(fileread(path));
            catch exception
                error('OpenOceanKraken:InvalidJSON', 'Cannot decode %s: %s', path, exception.message);
            end
            if ~isstruct(config)
                error('OpenOceanKraken:InvalidJSON', 'OOK JSON root must be an object.');
            end
            obj.RawConfig = config;
            obj.InputPath = path;
            obj.applyConfig(config);
        end

        function setBiologicalAttenuation(obj, layers, attenuationUnit)
            if nargin < 3
                attenuationUnit = 'dB/lambda';
            end
            allowedUnits = {
                'dB/m/kHz', 'Params_Lose', 'dB/m', ...
                'Nepers/m', 'Quality_Factor', 'dB/lambda'};
            if isstring(attenuationUnit) && isscalar(attenuationUnit)
                attenuationUnit = char(attenuationUnit);
            end
            if ~ischar(attenuationUnit) || ~isrow(attenuationUnit) || ...
                    ~any(strcmp(attenuationUnit, allowedUnits))
                error('OpenOceanKraken:InvalidAttenuationUnit', ...
                    'attenuationUnit must be one of the supported exact unit strings.');
            end
            if ~isnumeric(layers) || ~isreal(layers) || ...
                    ndims(layers) ~= 2 || size(layers, 2) ~= 5 || ...
                    size(layers, 1) > 200 || any(~isfinite(layers), 'all')
                error('OpenOceanKraken:InvalidBiologicalLayer', ...
                    'layers must be a finite real N-by-5 matrix with N <= 200.');
            end
            if any(layers(:, 1) > layers(:, 2)) || ...
                    any(layers(:, 3) <= 0) || ...
                    any(layers(:, 4) <= 0) || ...
                    any(layers(:, 5) < 0)
                error('OpenOceanKraken:InvalidBiologicalLayer', ...
                    'Require Z1 <= Z2, f0 > 0, Q > 0, and a0 >= 0.');
            end

            encodedLayers = cell(1, size(layers, 1));
            for index = 1:size(layers, 1)
                encodedLayers{index} = struct( ...
                    'Z1', layers(index, 1), ...
                    'Z2', layers(index, 2), ...
                    'f0', layers(index, 3), ...
                    'Q', layers(index, 4), ...
                    'a0', layers(index, 5));
            end
            obj.Attenuation = struct( ...
                'AttenuationUnit', attenuationUnit, ...
                'OceanAbsorptionModel', 'Biological', ...
                'BiologicalLayers', {encodedLayers});
        end

        function jsonPath = Write(obj, cachePath)
            cachePath = char(cachePath);
            [parent, ~, extension] = fileparts(cachePath);
            if strcmpi(extension, '.json')
                jsonPath = cachePath;
                if ~isempty(parent) && ~isfolder(parent)
                    mkdir(parent);
                end
            else
                if ~isfolder(cachePath)
                    mkdir(cachePath);
                end
                jsonPath = fullfile(cachePath, 'OpenOceanKraken_input.json');
            end
            config = obj.makeConfig();
            try
                text = jsonencode(config, 'PrettyPrint', true);
            catch
                text = jsonencode(config);
            end
            fid = fopen(jsonPath, 'w', 'n', 'UTF-8');
            if fid < 0
                error('OpenOceanKraken:WriteFailed', 'Cannot create JSON file: %s', jsonPath);
            end
            cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
            if fprintf(fid, '%s\n', text) < 0
                error('OpenOceanKraken:WriteFailed', 'Cannot write JSON file: %s', jsonPath);
            end
        end

        function resultRoot = run(obj, cachePath, varargin)
            parser = inputParser;
            addParameter(parser, 'Velocity', obj.VelocityEnabled, @(x) islogical(x) && isscalar(x));
            addParameter(parser, 'Mod', false, @(x) islogical(x) && isscalar(x));
            addParameter(parser, 'ModOnly', false, @(x) islogical(x) && isscalar(x));
            addParameter(parser, 'Json', false, @(x) islogical(x) && isscalar(x));
            parse(parser, varargin{:});
            options = parser.Results;
            if isempty(obj.InputPath) || ~isfile(obj.InputPath)
                error('OpenOceanKraken:MissingInput', 'Call loadEnv or loadJson before run.');
            end
            if ~isscalar(obj.NumThreads) || obj.NumThreads < 1 || obj.NumThreads ~= floor(obj.NumThreads)
                error('OpenOceanKraken:InvalidThreads', 'NumThreads must be a positive integer.');
            end
            cachePath = char(cachePath);
            if ~isfolder(cachePath)
                mkdir(cachePath);
            end
            resultRoot = fullfile(cachePath, 'OpenOceanKraken_result');
            command = sprintf('"%s" "%s" --output "%s" --threads %d', ...
                obj.Executable, obj.InputPath, resultRoot, obj.NumThreads);
            if options.Velocity
                command = [command ' --velocity']; %#ok<AGROW>
            end
            if options.Mod
                command = [command ' --mod']; %#ok<AGROW>
            end
            if options.ModOnly
                command = [command ' --mod-only']; %#ok<AGROW>
            end
            if options.Json
                command = [command ' --json']; %#ok<AGROW>
            end
            [status, output] = system(command);
            if status ~= 0
                error('OpenOceanKraken:RunFailed', ...
                    'OpenOceanKraken exited with status %d:\n%s', status, output);
            end
            if options.ModOnly
                obj.requireFile([resultRoot '.mod']);
            elseif options.Velocity
                obj.requireFile([resultRoot '_P.shd']);
                obj.requireFile([resultRoot '_V.shd']);
                obj.requireFile([resultRoot '_H.shd']);
            else
                obj.requireFile([resultRoot '.shd']);
            end
            if options.Mod
                obj.requireFile([resultRoot '.mod']);
            end
            if options.Json
                obj.requireFile([resultRoot '.json']);
            end
            obj.LastResultRoot = resultRoot;
        end

        function field = getPressure(obj, location)
            root = obj.resolveRoot(location);
            velocityPath = [root '_P.shd'];
            if isfile(velocityPath)
                field = read_shd(velocityPath);
            else
                field = read_shd([root '.shd']);
            end
        end

        function field = getVerticalVelocity(obj, location)
            field = read_shd([obj.resolveRoot(location) '_V.shd']);
        end

        function field = getHorizontalVelocity(obj, location)
            field = read_shd([obj.resolveRoot(location) '_H.shd']);
        end

        function modes = getModes(obj, location)
            modes = read_mod([obj.resolveRoot(location) '.mod']);
        end

        function axesHandle = plotPressure(obj, location, varargin)
            axesHandle = plotshd(obj.getPressure(location), varargin{:});
        end

        function figureHandle = plotModes(obj, location, varargin)
            figureHandle = plotmode(obj.getModes(location), varargin{:});
        end
    end

    methods (Access = private)
        function applyConfig(obj, config)
            if isfield(config, 'Title'), obj.Title = config.Title; end
            if isfield(config, 'freqinfo')
                if isfield(config.freqinfo, 'freq'), obj.Frequency = config.freqinfo.freq; end
                if isfield(config.freqinfo, 'freqvec')
                    obj.FrequencyVector = config.freqinfo.freqvec;
                end
            end
            if isfield(config, 'AttenUnit'), obj.Attenuation = config.AttenUnit; end
            if isfield(config, 'Pos')
                position = config.Pos;
                if isfield(position, 'GridType'), obj.GridType = position.GridType; end
                if isfield(position, 'SrcDepth'), obj.SourceDepth = position.SrcDepth; end
                if isfield(position, 'RecvDepth')
                    obj.ReceiverDepth = obj.expandGrid(position.RecvDepth, 'NRz');
                end
                if isfield(position, 'RecvRange')
                    obj.ReceiverRange = obj.expandGrid(position.RecvRange, 'NRr');
                end
            end
            if isfield(config, 'MLimit'), obj.MLimit = config.MLimit; end
            if isfield(config, 'RProf'), obj.ProfileRange = config.RProf; end
            if isfield(config, 'sspInput'), obj.SSP = config.sspInput; end
            if isfield(config, 'ReflectionCoef'), obj.ReflectionCoef = config.ReflectionCoef; end
            if isfield(config, 'SBP'), obj.SBP = config.SBP; end
            if isfield(config, 'is_Velocity'), obj.VelocityEnabled = config.is_Velocity; end
            if isfield(config, 'cLow'), obj.PhaseSpeed(1) = config.cLow; end
            if isfield(config, 'cHigh'), obj.PhaseSpeed(2) = config.cHigh; end
            if isfield(config, 'Rmax'), obj.Rmax = config.Rmax; end
            if isfield(config, 'SourceType'), obj.SourceType = config.SourceType; end
            if isfield(config, 'RunMode'), obj.RunMode = config.RunMode; end
            if isfield(config, 'CoherenceType'), obj.CoherenceType = config.CoherenceType; end
            if isfield(config, 'ModeType'), obj.ModeType = config.ModeType; end
        end

        function config = makeConfig(obj)
            if isempty(fieldnames(obj.RawConfig))
                config = struct();
            else
                config = obj.RawConfig;
            end
            config.Title = obj.Title;
            frequencyVector = obj.FrequencyVector;
            if isempty(frequencyVector)
                frequencyVector = obj.Frequency;
            end
            config.freqinfo = struct('Nfreq', numel(frequencyVector), 'freq', obj.Frequency);
            if numel(frequencyVector) > 1
                config.freqinfo.freqvec = frequencyVector;
            end
            config.AttenUnit = obj.Attenuation;
            if ~isfield(config, 'Pos'), config.Pos = struct(); end
            config.Pos.GridType = obj.GridType;
            config.Pos.SrcDepth = obj.forceNumericArray(obj.SourceDepth);
            config.Pos.RecvRange = obj.compactGrid(obj.ReceiverRange, 'NRr');
            config.Pos.RecvDepth = obj.compactGrid(obj.ReceiverDepth, 'NRz');
            if ~isfield(config.Pos, 'RecvAzim')
                config.Pos.RecvAzim = struct('start', 0, 'end', 0, 'NRo', 1);
            end
            config.MLimit = obj.MLimit;
            config.NProf = numel(obj.ProfileRange);
            config.RProf = obj.forceNumericArray(obj.ProfileRange);
            if ~isfield(config, 'hasModePos'), config.hasModePos = false; end
            config.sspInput = obj.forceSspArray(obj.SSP);
            config.ReflectionCoef = obj.ReflectionCoef;
            config.SBP = obj.SBP;
            config.is_Velocity = obj.VelocityEnabled;
            config.cLow = obj.PhaseSpeed(1);
            config.cHigh = obj.PhaseSpeed(2);
            if obj.Rmax > 0
                config.Rmax = obj.Rmax;
            elseif ~isempty(obj.ReceiverRange)
                config.Rmax = max(obj.ReceiverRange);
            else
                config.Rmax = 0;
            end
            config.SourceType = obj.SourceType;
            config.RunMode = obj.RunMode;
            config.CoherenceType = obj.CoherenceType;
            config.ModeType = obj.ModeType;
        end

        function root = resolveRoot(obj, location)
            if nargin < 2 || isempty(location)
                root = obj.LastResultRoot;
            elseif isfolder(location)
                root = fullfile(char(location), 'OpenOceanKraken_result');
            else
                root = char(location);
            end
            if isempty(root)
                error('OpenOceanKraken:MissingResult', 'No result root is available.');
            end
        end
    end

    methods (Static, Access = private)
        function result = expandGrid(value, countName)
            if isstruct(value) && isfield(value, 'start') && isfield(value, 'end') && isfield(value, countName)
                result = linspace(value.start, value.end, value.(countName));
            else
                result = value;
            end
        end

        function result = compactGrid(values, countName)
            values = values(:).';
            if isempty(values)
                result = [];
            else
                result = struct('start', values(1), 'end', values(end), countName, numel(values));
            end
        end

        function result = forceNumericArray(values)
            values = values(:).';
            if isscalar(values)
                result = {values};
            else
                result = values;
            end
        end

        function result = forceSspArray(value)
            if isempty(value)
                result = {};
                return;
            end
            if iscell(value)
                result = value;
                return;
            end
            result = cell(1, numel(value));
            for index = 1:numel(value)
                area = value(index);
                if isfield(area, 'layers') && isstruct(area.layers)
                    area.layers = num2cell(area.layers);
                end
                result{index} = area;
            end
        end

        function requireFile(path)
            if ~isfile(path)
                error('OpenOceanKraken:MissingResult', 'Expected output file was not created: %s', path);
            end
        end
    end
end
