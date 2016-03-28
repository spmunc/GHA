function [X_centered] = zero_mean(X)
% Centers X such that the columns have 0 mean

    X_centered = bsxfun(@minus,X,mean(X));

end

