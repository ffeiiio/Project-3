function dydt = net_sys_rhs(xx, par)
    % Wilson-Cowan 网络模型 (Network Wilson-Cowan Model)
    % xx: 包含当前变量和时滞变量的矩阵
    % par: 参数向量 [c1, c2, c3, c4, P, Q, beta, epsilon, tau, rho]

    % 解析参数
    c1 = par(1); c2 = par(2); c3 = par(3); c4 = par(4);
    P = par(5); Q = par(6); beta = par(7);
    epsilon = par(8); 
    tau = par(9);    % τ 代表 u(t - τ)
    rho = par(10);   % ρ 代表 u(t - ρ)

    % Wilson-Cowan 传递函数
    f = @(x) 1 ./ (1 + exp(-beta * x));

    % 计算系统动力学
    dydt = [...
        -xx(1,1) + f(c1 * xx(1,2) + c2 * xx(2,2) + epsilon * xx(1,3) + P);  % u 变量
        -xx(2,1) + f(c3 * xx(1,2) + c4 * xx(2,2) + Q)                       % v 变量
    ];
end
