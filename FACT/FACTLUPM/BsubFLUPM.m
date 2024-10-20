% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------                                                                                                               ------------------------------
% ------------------------------                                                BsubFLUPM                                                      ------------------------------
% ------------------------------                                                                                                               ------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



global L;
global P;
global Q;
global ansb;



function [x] = BsubFLUPM(A, b, L, P, Q)

  % L * y = b
  [n,n]=size(L);
  b = P*b;
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


  % U * z = y
  [n,n]=size(A);
  z = size(A,1);
  z = zeros(z,1);
  z(n)=y(n)/A(n,n);
  for i=n-1:-1:1
    S=0;
    for j=i+1:n
      S=S+A(i,j)*z(j);
    endfor
    z(i)=(y(i)-S)/A(i,i);
  endfor

  disp("Z:");
  disp(z);

  x = Q*z;
endfunction

