function [stacked] = Stack3D(X)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

[n,m,l] = size(X);

stacked = zeros(n*l,m);


for i=1:l
    inds = ((i-1)*n+1):(i*n);
    stacked(inds,:) = X(:,:,i);
end

