clear; clc; close all;
load Matrix;

[n] = size(Matrix, 1);
G = Matrix(1 : n - 2, 1 : n - 2);
b  = Matrix(1:n - 2, n - 1);

[x, iter] = jacobi(eye(n - 2) - G, b, 10000, 1e-6);

disp("Probabilitatile sunt:")
disp(x);

disp("Numarul necesar de iteratii:")
disp(iter)