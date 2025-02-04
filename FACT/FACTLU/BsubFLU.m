% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------                                                                                                               ------------------------------
% ------------------------------                                                 BsubFLU                                                       ------------------------------
% ------------------------------                                                                                                               ------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



global L;
global ansb;



function [x] = BsubFLU(A, b, L)

  % L * y = b
  [n,n]=size(L);
  y = zeros(n,1);
  y(1) = b(1) / L(1,1);
  for i=2:n
    S=0;
    for j=1:i-1
      S=S+L(i,j)*y(j);
    endfor
    y(i)=(b(i)-S)/L(i,i);
  endfor


  disp("Y:");
  disp(y);


  % U * x = y
  [n,n]=size(A);
  x(n)=y(n)/A(n,n);
  for i=n-1:-1:1
    S=0;
    for j=i+1:n
      S=S+A(i,j)*x(j);
    endfor
    x(i)=(y(i)-S)/A(i,i);
  endfor
endfunction

