function numElements = fnmemb(fileName, elementSize)
%% FNMEMB  Calculate number of elements in a binary file.
%
% Syntax:
%   numElements = fnmemb(fileName)
%   numElements = fnmemb(fileName, elementSize)
%
% Inputs:
%   fileName     - Path to the binary file
%   elementSize  - Size of each element in bytes (default is 4, i.e., float32)
%
% Output:
%   numElements  - Number of elements in the file

% Default element size if not provided
if nargin < 2
    elementSize = 4;  % float32 default
end

% Validate elementSize
if elementSize <= 0
    error('Element size must be a positive integer.');
end

% Get total file size in bytes
fileSizeInBytes = fsize(fileName);

% Calculate number of elements
numElements = fileSizeInBytes / elementSize;

% Ensure result is an integer
if mod(numElements, 1) ~= 0
    warning('File size is not a multiple of element size. Result is not an integer.');
end

end
