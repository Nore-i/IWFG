% y = [2, 8, 6, 1;
%     5, 4, 9, 1;
%     7, 2, 3, 8;
%     9, 6, 1, 7;
%     3, 5, 7, 4;
%     5, 5, 5, 5];
% y = [2, 8, 6;
%     5, 4, 9;
%     7, 2, 3;
%     9, 6, 1;
%     3, 5, 7;
%     5, 5, 5];
y=[ 0.1085,    3.7468;
    0.3385,    3.6288;
    0.0643,    4.0344;
    0.9837,    2.7669;
    0.3385,    3.0072;
    0.4475,    2.8060;
    1.0000,    2.6218;
    0.3385,    3.2356];

% y = 0.1 .* y;
ref = [1.2000, 4.8413]; % WFG convention: maximization problems
% y = -y; % Our convention: minimization problems, so ref > y. 

result = stk_dominatedhv(y, ref, 0)
% indices = stk_IWFG(y, ref, 4)
% stk_IWFG_mex (y, 2)