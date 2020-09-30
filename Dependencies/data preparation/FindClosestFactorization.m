function factors = FindClosestFactorization(y)
%FINDCLOSEFACTORS Given y, finds x1 and x2 such that x1*x2 = y and |x1-x2| is minimized.
%   If y is not an integer, it will be rounded towards the nearest one.
%
%   Examples:
%       Input: 12   Output: [3 4]
%       Input: 39   Output: [3 13]
%       Input: 4    Output: [2 2]
%       Input: 7    Output: [1 7]
%
%   Author: Luke Gane
%	Version: 1.0
%   Last updated: 2016-05-21

y = round(y(:));

firstFactor = arrayfun(@FindFirstFactor, y);

secondFactor = y./firstFactor;

factors = [firstFactor secondFactor];

end

function firstFactor = FindFirstFactor(value)

firstFactor = floor(sqrt(value));

while (mod(value, firstFactor) ~= 0)
    firstFactor = firstFactor - 1;
end

end