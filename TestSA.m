clear ; close all;
A = [2,1,1,3;1,3,1,2];
b = [5;3];
c = [6;8;5;9];
[b_bar, eps_bar, A_bar, c_bar, B, NonB] = SimplexAlg(A, b, c)