## Copyright Horia Mercan 2023

## @deftypefn {} {@var{x, iter} =} power_method (@var{A}, @var{b}, @var{max_iter}, @var{tol})
##
## @seealso{}
## @end deftypefn

function [x, iter] = jacobi(A, b, max_iter, tol)
  [n] = size(A, 1);
  D = diag(diag(A));
  P = D - A;
  
  disp("The spectral radius of the iteration matrix")
  max(abs(eig(inv(D) * P)))
  G = inv(D) * P;
  b = inv(D) * b;
  x0 = zeros(n, 1);
  
  for iteration = 1 : max_iter
    iter = iteration;
    
    x = G * x0 + b;
    if (norm(x - x0) < tol)
      break;
    endif
    
    x0 = x;
  endfor
  
  decimals = round(-log(tol)/log(10));
  x = round(x * 10^decimals);
  x = x / 10^decimals;
endfunction