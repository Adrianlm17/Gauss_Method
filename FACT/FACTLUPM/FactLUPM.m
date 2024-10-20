% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------                                                                                                               ------------------------------
% ------------------------------                                                 FactLUPM                                                      ------------------------------
% ------------------------------                                                                                                               ------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



global L;
global Q;
global P;

function [A, L, P, Q] = FactLUPM(A)
  % Inicializar L, P, Q como matrices identidad
  global L = eye(size(A));
  L = eye(size(A));
  global P = eye(size(A));
  P = eye(size(A));
  global Q = eye(size(A));
  Q = eye(size(A));

  [m, n] = size(A);

  % Vector de índice de columnas
  w = 1:n;

  % Eliminar variables para triangular la matriz
  for k = 1:n-1
    % Encontrar el máximo pivote en la submatriz A(k:m, k:n)
    [pivot_value, row_pivot] = max(abs(A(k:m, k:n)(:)));
    [row_offset, col_offset] = ind2sub([m-k+1, n-k+1], row_pivot);
    row_pivot = row_offset + k - 1;
    col_pivot = col_offset + k - 1;

    % Intercambiar filas si es necesario
    if row_pivot ~= k
      A([k, row_pivot], :) = A([row_pivot, k], :);
      P([k, row_pivot], :) = P([row_pivot, k], :); % Actualizar matriz de permutación de filas
    endif

    % Intercambiar columnas si es necesario
    if col_pivot ~= k
      A(:, [k, col_pivot]) = A(:, [col_pivot, k]);
      Q(:, [k, col_pivot]) = Q(:, [col_pivot, k]); % Actualizar matriz de permutación de columnas
      w([k, col_pivot]) = w([col_pivot, k]); % Actualizar vector de índices
    endif

    % Eliminar elementos debajo del pivote y actualizar la matriz L
    for i = k+1:m
      L(i, k) = A(i, k) / A(k, k); % Guardar el multiplicador en L
      A(i, k:n) = A(i, k:n) - L(i, k) * A(k, k:n); % Restar la fila escalada
    endfor
  endfor
endfunction

