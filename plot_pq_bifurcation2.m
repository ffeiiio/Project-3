clear; clc; close all;

% 添加 DDE-Biftool 相关路径
addpath(genpath('路径到ddebiftool/'));

% 设置 Wilson-Cowan 参数
beta = 1;
c1 = 10; 
c2 = -10; 
c3 = 10; 
c4 = 2;

% 设定 P 和 Q 的范围（点数适当增加）
P_values = linspace(-10, 10, 150);  
Q_values = linspace(-10, 10, 150);

% 存储 Hopf 分岔点 和 Saddle-Node 分岔点
hopf_P = []; hopf_Q = [];
saddle_P = []; saddle_Q = [];

% 遍历 P, Q 计算分岔
for P = P_values
    for Q = Q_values
        % 计算稳态
        y0 = [0.5; 3]; % 初始猜测值
        steady_state_func = @(y) wilson_cowan_ode(0, y, P, Q, beta, c1, c2, c3, c4);
        y_star = fsolve(steady_state_func, y0, optimoptions('fsolve', 'Display', 'off'));

        % 计算雅可比矩阵
        J = [
            -1 + c1 * beta * y_star(1) * (1 - y_star(1)), c2 * beta * y_star(1) * (1 - y_star(1));
            c3 * beta * y_star(2) * (1 - y_star(2)), -1 + c4 * beta * y_star(2) * (1 - y_star(2))
        ];

        % 计算特征值
        eigenvalues = eig(J);

        % **Hopf 分岔条件**
        % - Tr(J) ≈ 0 （迹接近 0）
        % - det(J) > 0 （行列式大于 0）
        % - 存在一对复特征值，实部接近 0
        if abs(trace(J)) < 0.1 && det(J) > 0.1
            hopf_P = [hopf_P, P];
            hopf_Q = [hopf_Q, Q];
        end

        % **Saddle-Node 分岔条件**
        % - det(J) ≈ 0 （行列式接近 0）
        if abs(det(J)) < 0.1
            saddle_P = [saddle_P, P];
            saddle_Q = [saddle_Q, Q];
        end
    end
end

% **绘制 Hopf 和 Saddle-Node 分岔图**
figure; hold on;
scatter(hopf_P, hopf_Q, 'r', 'filled');  % Hopf 分岔点（红色）
scatter(saddle_P, saddle_Q, 'b', 'filled');  % Saddle-Node 分岔点（蓝色）

xlabel('P');
ylabel('Q');
title('Bifurcation Diagram in P-Q Plane');
grid on;

legend({'Hopf Bifurcation (Red)', 'Saddle-Node Bifurcation (Blue)'}, 'Location', 'NorthEast');