n1 =1200;
n2 = 900;
s=300;
r = 2;
rng("default");
SU = generateSketchingMatrix(n1,s);
SV = generateSketchingMatrix(n2,s);


A= diag(rand(1,n1));
B= diag(rand(1,n2));
C1 = rand(n1,r);
C2 = rand(n2,r);
 
C = C1*C2';

X_lyap = lyap(A,B,-C);





k=1;
maxit =100;
tol = 1e-8;
p=1;
[X1,X2] = STArnoldi_computation(A, B, C1, C2, SU, SV, k, maxit, tol, p);

Xd = X1*X2';

error2 = norm(Xd-X_lyap, 'fro')/norm(X_lyap,'fro');
errorr = norm(A*Xd+Xd*B-C1*C2','fro')/norm(C1*C2','fro');
disp(['Approximation error X_d - X: ', num2str(error2)]);
disp(['residue: ', num2str(errorr)]);
















function S = generateSketchingMatrix(n, s)
    % n: dimension of the original matrix (n x n)
    % s: number of rows in the sketching matrix (s < n)
    rng("default")

    % Generate E: Diagonal matrix with Rademacher random variables (-1 or 1)
    E = diag(randsrc(n, 1, [-1, 1]));

    % N: Discrete cosine transform matrix
    N = dct(eye(n));

    % Generate D: Matrix to randomly select s rows
    rowIndices = randperm(n, s); % Randomly select s unique indices from 1 to n
    D = zeros(s, n);
    for i = 1:s
        D(i, rowIndices(i)) = 1;
    end

    % Calculate the sketching matrix S
    S = sqrt(s / n) * D * N * E;
end

