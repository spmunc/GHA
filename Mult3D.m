function [mult] = Mult3D(X,Y)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
[xx,xy,xz] = size(X);
[yx,yy,yz] = size(Y);
assert(xy == yx, 'Inner Dimensions Must Agree')
assert(xz == yz, '3rd Dimensional Length Must Agree');
mult = zeros(xx,yy,yz);
for i = 1:yz
    mult(:,:,i) = X(:,:,i)*Y(:,:,i);
end

end

