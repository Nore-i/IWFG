y = [2, 8, 6, 1;
    5, 4, 9, 1;
    7, 2, 3, 8;
    9, 6, 1, 7;
    3, 5, 7, 4;
    5, 5, 5, 5];
% y = [2, 8, 6;
%     5, 4, 9;
%     7, 2, 3;
%     9, 6, 1;
%     3, 5, 7;
%     5, 5, 5];

y = 0.1 .* y;
ref = [0, 0, 0, 0]; % WFG convention: maximization problems
y = -y; % Our convention: minimization problems, so ref > y. 

result = stk_dominatedhv(y, ref, 0);
indices = stk_IWFG(y, ref, 4)