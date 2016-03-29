function [diags] = Diag3D(X)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[~,y,z] = size(X);
diags = zeros(z,y);

for i=1:z
    diags(i,:) = diag(X(:,:,i))';
end

end

