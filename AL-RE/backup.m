
C = rand(100, 50); % An example matrix
r = 10;            % Desired rank of the approximation

[C1, C2] = LowRankApproximation(C, r);

% Verify the approximation
C_approx = C1 * C2';
error = norm(C - C_approx, 'fro');
disp(['Approximation error (Frobenius norm): ', num2str(error)]);
function [C1, C2] = LowRankApproximation(C, r)
    % C: The original matrix to be approximated
    % r: The desired rank of the approximation
    
    % Compute the SVD of C
    [U, S, V] = svd(C, 'econ');
    
    % Select the first r singular values/vectors
    Ur = U(:, 1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:, 1:r);
    
    % Compute the square root of the singular values matrix
    SqrtSr = sqrt(Sr);
    
    % Define C1 and C2 for the low-rank approximation C ¡Ö C1*C2'
    C1 = Ur * SqrtSr;
    C2 = Vr * SqrtSr;
end