function [ numLines ] = GetNumberOfLines( filename )
    fid = fopen(filename, 'rb');
    %# Get file size.
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    frewind(fid);
    %# Read the whole file.
    data = fread(fid, fileSize, 'uint8');
    %# Count number of line-feeds and increase by one.
    numLines = sum(data == 10) + 1;
    fclose(fid);
end
