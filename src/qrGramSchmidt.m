%Calculo de las matrices Q y R utilizando Gram-Schmidt a partir de la matriz A.
function [Q, R] = qrGramSchmidt(A)
[~, cols] = size(A);
for x = 1:cols
   col = A(:, x);
   for y = 1:(x - 1)
        Qt = Q(:, y);
        R(y, x) = Qt' * A(:, x);
        col -= R(y, x) * Qt;
   end
   R(x, x) = norm(col);
   Q(:, x) = col / R(x, x);
end
