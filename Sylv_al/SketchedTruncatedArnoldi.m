function [X1, X2] = SketchedTruncatedArnoldi(A, B, C1, C2, SU, SV, k, maxit, tol, p)
    % Initialize variables
    [U1, ~] = qr(C1, 0);
    [V1, ~] = qr(C2, 0);
    [QU1, beta1] = qr(SU * C1, 0);
    [QV1, beta2] = qr(SV * C2, 0);
    TU1 = beta1;
    TV1 = beta2;
    Hhat = zeros(size(A, 1), size(B, 1));
    Ghat = zeros(size(A, 1), size(B, 1));

    for d = 1:maxit
        % Compute U and V
        U = A * U1(:, d);
        V = B' * V1(:, d);

        % Orthogonalization
        for i = max(1, d - k + 1):d
            U = U - U1(:, i) * (U1(:, i)' * U);
            V = V - V1(:, i) * (V1(:, i)' * V);
        end

        % Compute skinny QRs
        [U1(:, d+1), h] = qr(U, 0);
        [V1(:, d+1), g] = qr(V, 0);

        % Update QRs
        [QU1, ~] = qr(SU * [U1(:, 1:d), U1(:, d+1)], 0);
        [QV1, ~] = qr(SV * [V1(:, 1:d), V1(:, d+1)], 0);

        % Update Hhat and Ghat
        tauD1Inv = inv(TU1(:, end)); % Inverse of the last element of TU1
        Hnew = (Hhat + h * beta1' * tauD1Inv) * tauD1Inv + TU1 * H * tauD1Inv + h * h(d+1, d) * tauD1Inv;
        Hhat = [Hhat + h * beta1' * tauD1Inv, Hnew; ...
                tauD1Inv * h * beta1' * tauD1Inv, tauD1Inv * (h * beta1' * TU1 + h(d+1, d+1) * tauD1Inv)];
        % Similar update for Ghat

        if mod(d, p) == 0
            % Solve the reduced Sylvester equation
            Y = lyap(Hhat, Ghat, QU1(:, 1:d)' * C1 * C2' * QV1(:, 1:d));
            rho = sqrt(norm(h * Y, 'fro')^2 + norm(Y * g', 'fro')^2);

            if rho < tol
                break;
            end
        end
    end

    % Compute low-rank factors
    [Y1, Y2] =
    % ...continuing from the previous code

    % Implement LowRankFactors function to compute low-rank factors of Y
    [Y1, Y2] = LowRankFactors(Y);

    % Retrieve solutions
    X1 = U1 * inv(TU1) * Y1;
    X2 = V1 * inv(TV1) * Y2;
end


