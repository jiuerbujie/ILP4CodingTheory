function [b_bar, eps_bar, A_bar, c_bar, B, NonB] = SimplexAlg(A, b, c)
% Solve max c'x s.t. Ax <= b,x>0  by Simplex Algorithm
% Output:
% b_bar: The value of basis elements.
% eps_bar: The optimal value of object function
% A_bar: The A_bar in dictionary.
% B: basis index
% N: non-basis index

% Check feasibility
if sum(b>0) ~=length(b)
    disp(['Simplex Algorithm infeasible']);
end
c = [c;zeros(size(A,1),1)];
% Initial dictionary
A_bar = A;
b_bar = b;
NonB = 1:size(A,2);
B = 1+size(A,2):size(A,1)+size(A,2);
c_B = c(B);
c_N = c(NonB);
eps_bar = c_B' * b_bar;
c_bar = (c_B'*A_bar - c_N')';


% Main Iteration
while AlgStopCheck(c_bar, b_bar)
% Find pivot element
    [~,j] = min(c_bar);
    idx_tmp = find(A_bar(:,j) > 0);
    b_div_a_tmp = min(b_bar(idx_tmp)./A_bar(idx_tmp,j));
    i = find(b_bar./A_bar(:,j) == b_div_a_tmp);
    %i = find(b_bar == min(b_bar(idx_tmp)./A_bar(idx_tmp,j)));
    A_back = A_bar;
    pivot = A_bar(i,j);
    tmp = NonB(j);
    NonB(j) = B(i);
    B(i) = tmp;
    
% Elimilation
    LargeA = [A_bar,b_bar;c_bar',eps_bar];  % The whole dictionary
    LargeA_back = LargeA;
    NonPivotRow = 1:length(B)+1;
    NonPivotRow(i) = [];
    for k = NonPivotRow
        LargeA(k,:) = LargeA(k,:)-LargeA(k,j)/pivot*LargeA(i,:);
    end
    LargeA(i,:) = LargeA_back(i,:) / pivot;
    LargeA(:,j) = LargeA_back(:,j) / pivot*-1;
    LargeA(i,j) = LargeA_back(i,j) / pivot / pivot;
    
    A_bar = LargeA(1:length(B),1:length(NonB));
    b_bar = LargeA(1:length(B),end);
    c_bar = LargeA(end,1:length(NonB))';
    eps_bar = LargeA(end,end);
%     
%     A_bar(i,:) = A_back(i,:) / pivot;
%     A_bar(:,j) = A_back(:,j) / pivot*-1;
%     A_bar(i,j) = A_back(i,j) / pivot / pivot;
    
end
return

function r = AlgStopCheck(c_bar, b_bar)
r = 1;
if sum(c_bar > 0 ) ==length(c_bar)
    disp('Optimal vaule reached!');
    r = 0;
    return
end
if sum(b_bar > 0 )~= length(b_bar)
    disp('LP becomes infeasible');
    r = 0;
    return
end
return