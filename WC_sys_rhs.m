function dydt = WC_sys_rhs(xx, par)
    % Wilson-Cowan 模型的右端项
    % xx: 包含当前变量和时滞变量的矩阵
    % par: 参数向量 [c1, c2, c3, c4, P, Q, beta, tau1, tau2]

    dydt = [...
        -xx(1,1) + 1./(1 + exp(-par(7) * (par(1) * xx(1,2) + par(2) * xx(2,3) + par(5))));
        -xx(2,1) + 1./(1 + exp(-par(7) * (par(3) * xx(1,3) + par(4) * xx(2,2) + par(6))))
    ];
end