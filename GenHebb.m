function [W, Y, iters, innerWs, checks, errors] = GenHebb(X, mu, maxiter, tol, chkpnts)
% GenHebb runs Sanger's Generalized Hebbian Algorithm to learn all 
% principle eigenvectors of the correlation matrix of X. 
%   Inputs:
%       X: mxn matrix where the inputs are rows
%       mu: learning rate
%       maxiter: maximum number of learning steps
%       tol: if ||W*W' - I|| < tol (matrix 2-norm), then the learning stops
%       chkpnts: how many times to check in on and record progress
%                     
%   Outputs: 
%       W: weight matrix, where rows are the learned principle eigenvectors
%       Y: output of network, which are principle components 
%       iters: number of learning steps actually performed
%       innerWs: values of W at checkpoints ( size: nxnxlength(checks)+2)
%       checks: learning steps at which innerWs are
%       errors: matrix 2 norm of W*W' - I

% Get size of X
[m,n] = size(X);

% Weight matrix is nxn
W_t = 2*rand(n,n) - 1;

% Initialize matrix to fill at checkpoints
innerWs = zeros(n,n,chkpnts+2);
innerWs(:,:,1) = W_t; % first entry is initial weight matrix
errors = zeros(chkpnts+2,1);
errors(1) = norm(W_t*W_t'-eye(n));

% checks
interval = floor(maxiter/chkpnts);
checks = interval:interval:maxiter;
checks = [0, checks, maxiter];
curr_chk = 2;

% Run the GHA algorithm
W_tplus1 = W_t;
for iter = 1:maxiter
    
    % Calculate current output
    Y_t = W_t*X'; % TODO: don't need all of this
    
    % Get the input and corresponding output for this step
    % t = randi([1,m],1,1); % pick the input
    t = mod(iter,m) + 1;
    x_t = X(t,:); % for consistency with formula
    y_t = Y_t(:,t); %for consistency with formula
    
    % update all weights
    for i = 1:n
        for j = 1:n
            summ = 0;
            for k=1:i
                summ = summ + W_t(k,j)*y_t(k);
            end
            W_tplus1(i,j) = W_t(i,j) + mu*y_t(i)*(x_t(j) - summ);
        end
    end
    
    % Prepare for next iteration
    W_t = W_tplus1;
    
    % If we are at a checkpoint, record W. 
    % Return if ||W*W' - I|| < tol
    if any(iter == checks)
        innerWs(:,:,curr_chk) = W_t;
        err = norm(W_t*W_t'-eye(n));
        errors(curr_chk) = norm(W_t*W_t'-eye(n));
        if err < tol
            % Rename for clarity
            W = W_tplus1;
            Y = Y_t;
            % fill end final W, errors, and iters 
            iters = iter;
            checks(end) = iter;
            errors(end) = norm(W_t*W_t'-eye(n));
            innerWs(:,:,end) = W_t;
            return
        end
        curr_chk = curr_chk + 1;
    end
    
end

% fill end final W and errors
innerWs(:,:,end) = W_t;
errors(end) = norm(W_t*W_t'-eye(n));

% Rename for clarity
W = W_tplus1;
Y = Y_t;
iters = maxiter;

end

