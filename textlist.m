
% Initialize a 5x5 cell array, each cell to hold a 2x2 matrix
blockMatrix = cell(5, 5);

% Fill each block with a 2x2 matrix (example: random matrices)
for i = 1:5
    for j = 1:5
        blockMatrix{i, j} = rand(2, 2); % Replace rand(2, 2) as needed
    end
end
% Assuming 'blockMatrix' is your large matrix composed of smaller blocks
% For example, a 10x10 matrix composed of 2x2 blocks

% First, ensure that 'blockMatrix' is square and invertible
if isequal(size(blockMatrix), [10 10]) && det(blockMatrix) ~= 0
    inverseMatrix = inv(blockMatrix);
else
    error('The entire block matrix is not square or invertible');
end
