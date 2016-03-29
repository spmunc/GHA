function [trans] = Trans3D(X)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

[x,y,z] = size(X);
trans = zeros(x,y,z);
assert(x == y, 'Matrices must be square');
for i = 1:z
    trans(:,:,i) = X(:,:,i)';
end

