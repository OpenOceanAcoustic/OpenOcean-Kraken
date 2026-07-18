function modes = read_mod(filename)
%READ_MOD Read every profile in an OpenOceanKraken MOD file.

errorId = 'OpenOceanKraken:InvalidMOD';
fid = fopen(filename, 'rb', 'ieee-le');
if fid < 0
    error(errorId, 'Cannot open MOD file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
info = dir(filename);
if isempty(info) || info.bytes < 4
    error(errorId, 'MOD file is too short.');
end

recl = readAt(fid, 0, 1, 'int32=>double', errorId);
if recl <= 0 || recl ~= floor(recl)
    error(errorId, 'Invalid MOD record length.');
end
recordBytes = 4 * recl;
if recordBytes > info.bytes || mod(info.bytes, recordBytes) ~= 0
    error(errorId, 'MOD fixed records are truncated.');
end

profiles = struct('profile_index', {}, 'profile_range', {}, 'depth', {}, ...
    'wavenumbers', {}, 'group_velocity', {}, 'mode_shapes', {}, ...
    'media_mesh_counts', {}, 'media_materials', {}, 'halfspaces', {});
recordIndex = 0;
title = '';
frequency = [];
while recordIndex * recordBytes < info.bytes
    base = recordIndex * recordBytes;
    embeddedRecl = readAt(fid, base, 1, 'int32=>double', errorId);
    if embeddedRecl ~= recl
        error(errorId, 'Inconsistent MOD record length in profile header.');
    end
    profileTitle = strtrim(deblank(char(readAt(fid, base + 4, 80, '*char', errorId).')));
    dims = readAt(fid, base + 84, 4, 'int32=>double', errorId);
    nfreq = dims(1);
    nmedia = dims(2);
    depthCount = dims(3);
    materialCount = dims(4); %#ok<NASGU>
    if nfreq <= 0 || nmedia <= 0 || depthCount <= 0 || any(dims ~= floor(dims))
        error(errorId, 'Invalid MOD profile dimensions.');
    end

    meshCounts = zeros(1, nmedia);
    materials = strings(1, nmedia);
    mediaBase = (recordIndex + 1) * recordBytes;
    for medium = 1:nmedia
        offset = mediaBase + (medium - 1) * 12;
        meshCounts(medium) = readAt(fid, offset, 1, 'int32=>double', errorId);
        materials(medium) = string(strtrim(deblank(char( ...
            readAt(fid, offset + 4, 8, '*char', errorId).'))));
    end

    frequencies = readAt(fid, (recordIndex + 3) * recordBytes, ...
        nfreq, 'double=>double', errorId);
    depth = readAt(fid, (recordIndex + 4) * recordBytes, ...
        depthCount, 'single=>double', errorId);
    modeRecord = recordIndex + 5;
    modeCount = readAt(fid, modeRecord * recordBytes, 1, 'int32=>double', errorId);
    if modeCount < 0 || modeCount ~= floor(modeCount)
        error(errorId, 'Invalid MOD mode count.');
    end
    halfspaces = readHalfspaces(fid, (modeRecord + 1) * recordBytes, errorId);

    modeShapes = complex(zeros(modeCount, depthCount, 'single'));
    for modeIndex = 1:modeCount
        raw = readAt(fid, (modeRecord + 1 + modeIndex) * recordBytes, ...
            2 * depthCount, 'single=>single', errorId);
        modeShapes(modeIndex, :) = reshape(complex(raw(1:2:end), raw(2:2:end)), 1, []);
    end

    modesPerRecord = max(1, floor(recl / 2));
    kRecordCount = ceil(modeCount / modesPerRecord);
    wavenumbers = complex(zeros(modeCount, 1, 'single'));
    first = 1;
    for kRecord = 1:kRecordCount
        take = min(modesPerRecord, modeCount - first + 1);
        raw = readAt(fid, (modeRecord + 1 + modeCount + kRecord) * recordBytes, ...
            2 * take, 'single=>single', errorId);
        wavenumbers(first:first + take - 1) = complex(raw(1:2:end), raw(2:2:end));
        first = first + take;
    end

    currentFrequency = frequencies(1);
    if ~isempty(frequency) && abs(currentFrequency - frequency) > ...
            eps(max(abs([currentFrequency, frequency]))) * 8
        error(errorId, 'MOD profiles have inconsistent frequencies.');
    end
    frequency = currentFrequency;
    if isempty(title)
        title = profileTitle;
    end
    profiles(end + 1) = struct( ... %#ok<AGROW>
        'profile_index', numel(profiles), ...
        'profile_range', NaN, ...
        'depth', depth, ...
        'wavenumbers', wavenumbers, ...
        'group_velocity', zeros(0, 1), ...
        'mode_shapes', modeShapes, ...
        'media_mesh_counts', meshCounts, ...
        'media_materials', materials, ...
        'halfspaces', halfspaces);
    if modeCount > 0
        extraKRecords = floor((2 * modeCount - 1) / recl);
    else
        extraKRecords = 0;
    end
    recordIndex = recordIndex + 8 + modeCount + extraKRecords;
end

if isempty(profiles)
    error(errorId, 'MOD file contains no profiles.');
end
modes = struct('title', title, 'frequency', frequency, ...
    'profiles', profiles, 'path', char(filename));
end

function halfspaces = readHalfspaces(fid, offset, errorId)
top = readHalfspace(fid, offset, errorId);
bottom = readHalfspace(fid, offset + 25, errorId);
halfspaces = struct('top', top, 'bottom', bottom);
end

function value = readHalfspace(fid, offset, errorId)
boundary = char(readAt(fid, offset, 1, '*char', errorId));
cpRaw = readAt(fid, offset + 1, 2, 'single=>double', errorId);
csRaw = readAt(fid, offset + 9, 2, 'single=>double', errorId);
densityDepth = readAt(fid, offset + 17, 2, 'single=>double', errorId);
value = struct('boundary', strtrim(boundary), ...
    'compressional_speed', complex(cpRaw(1), cpRaw(2)), ...
    'shear_speed', complex(csRaw(1), csRaw(2)), ...
    'density', densityDepth(1), 'depth', densityDepth(2));
end

function values = readAt(fid, offset, count, precision, errorId)
if fseek(fid, offset, 'bof') ~= 0
    error(errorId, 'Cannot seek to MOD record.');
end
[values, actual] = fread(fid, count, precision);
if actual ~= count
    error(errorId, 'Unexpected end of MOD file.');
end
end
