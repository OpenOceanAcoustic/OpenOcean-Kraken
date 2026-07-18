function field = read_shd(filename)
%READ_SHD Read an OpenOceanKraken SHD file.
% Values are returned in [source, range, depth] order.

errorId = 'OpenOceanKraken:InvalidSHD';
fid = fopen(filename, 'rb', 'ieee-le');
if fid < 0
    error(errorId, 'Cannot open SHD file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
info = dir(filename);
if isempty(info) || info.bytes < 4
    error(errorId, 'SHD file is too short.');
end

recl = readExact(fid, 1, 'int32=>double', errorId);
if recl <= 0
    error(errorId, 'Invalid SHD record length.');
end
recordBytes = 4 * recl;
if recordBytes > info.bytes || info.bytes < 10 * recordBytes
    error(errorId, 'SHD header records are truncated.');
end

seekExact(fid, 4, errorId);
title = strtrim(deblank(char(readExact(fid, 80, '*char', errorId).')));
seekExact(fid, recordBytes, errorId);
plotType = lower(strtrim(char(readExact(fid, 10, '*char', errorId).')));

seekExact(fid, 2 * recordBytes, errorId);
dims = readExact(fid, 7, 'int32=>double', errorId);
if any(dims <= 0) || any(dims ~= floor(dims))
    error(errorId, 'SHD dimensions are invalid.');
end
nfreq = dims(1);
nsz = dims(5);
nrz = dims(6);
nrr = dims(7);
frequencyHeader = readExact(fid, 1, 'single=>double', errorId);
attenuation = readExact(fid, 1, 'single=>double', errorId);

seekExact(fid, 3 * recordBytes, errorId);
frequencyVector = readExact(fid, nfreq, 'double=>double', errorId);
if isempty(frequencyVector)
    frequency = frequencyHeader;
else
    frequency = frequencyVector(1);
end
seekExact(fid, 7 * recordBytes, errorId);
sourceDepths = readExact(fid, nsz, 'single=>double', errorId);
seekExact(fid, 8 * recordBytes, errorId);
receiverDepths = readExact(fid, nrz, 'single=>double', errorId);
seekExact(fid, 9 * recordBytes, errorId);
receiverRanges = readExact(fid, nrr, 'single=>double', errorId);

lastRecordIndex = floor((info.bytes - 1) / recordBytes);
dataRecordCount = lastRecordIndex - 9;
if dataRecordCount <= 0 || mod(dataRecordCount, nsz) ~= 0
    error(errorId, 'SHD data records do not match the source count.');
end
depthCount = dataRecordCount / nsz;
isIrregular = startsWith(plotType, 'irregular');
if (isIrregular && depthCount ~= 1) || (~isIrregular && depthCount ~= nrz)
    error(errorId, 'SHD depth-record count is inconsistent with the grid type.');
end

values = complex(zeros(nsz, nrr, depthCount, 'single'));
for isz = 1:nsz
    for irz = 1:depthCount
        recordIndex = 10 + (isz - 1) * depthCount + (irz - 1);
        seekExact(fid, recordIndex * recordBytes, errorId);
        raw = readExact(fid, 2 * nrr, 'single=>single', errorId);
        row = complex(raw(1:2:end), raw(2:2:end));
        values(isz, :, irz) = reshape(row, 1, []);
    end
end

field = struct( ...
    'title', title, ...
    'plot_type', plotType, ...
    'frequency', frequency, ...
    'frequency_vector', frequencyVector, ...
    'attenuation', attenuation, ...
    'source_depths', sourceDepths, ...
    'receiver_ranges', receiverRanges, ...
    'receiver_depths', receiverDepths, ...
    'values', values, ...
    'path', char(filename));
end

function values = readExact(fid, count, precision, errorId)
[values, actual] = fread(fid, count, precision);
if actual ~= count
    error(errorId, 'Unexpected end of SHD file.');
end
end

function seekExact(fid, offset, errorId)
if fseek(fid, offset, 'bof') ~= 0
    error(errorId, 'Cannot seek to SHD record.');
end
end
