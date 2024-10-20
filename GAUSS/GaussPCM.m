% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------                                                                                                               ------------------------------
% ------------------------------                                          GAUSS Pivoting maximal                                               ------------------------------
% ------------------------------                                                                                                               ------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



function [A,b] = GaussPCM(A,b)
  [m, n] = size(A);

  if m ~= n
    error('Matrix A must be square!');
  endif

  % Index of columns
  w = 1:n;

  % Eliminate variables to triangulate the matrix
  for k = 1:n-1
    % Find the maximum pivot in the submatrix A(k:m, k:n)
    [max_val, row_pivot] = max(abs(A(k:m, k)));
    row_pivot = row_pivot + k - 1;

    if abs(max_val) == 0
      error('Matrix A is singular!');
    endif

    % Swap rows (if pivot is not on row k)
    if row_pivot ~= k
      A([k, row_pivot], :) = A([row_pivot, k], :);
      b([k, row_pivot], 1) = b([row_pivot, k], 1);
    endif

    % Find the maximum pivot in the remaining submatrix (swap columns)
    [~, col_pivot] = max(abs(A(k, k:n)));
    col_pivot = col_pivot + k - 1;

    if col_pivot ~= k
      A(:, [k, col_pivot]) = A(:, [col_pivot, k]);
      w([k, col_pivot]) = w([col_pivot, k]);
    endif

    % Remove terms below the pivot
    for i = k+1:m
      mult = A(i,k) / A(k,k);
      A(i,k:n) = A(i,k:n) - mult * A(k,k:n);
      b(i,1) = b(i,1) - mult * b(k,1);
    endfor
  endfor

  % Display vector w after column swaps
  disp('Vector w after column swaps:');
  disp(w);

  % Backward substitution resolution
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = (b(i,1) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
  endfor

  % Reorder x according to vector w
  x = x(w);
endfunction

