function d = WC_sys_mfderi(xx, par, nx, np, v)
    % 多阶导数计算 (mfderi) for Wilson-Cowan 模型
    % 输入:
    %   xx  - 变量历史状态 (dim × (ntau+1))
    %   par - 参数
    %   nx  - x 的偏导数阶数 (1, 2, 3)
    %   np  - 参数的偏导数阶数
    %   v   - 方向导数 (仅 nx=3 时使用)
    % 输出:
    %   d   - 高阶导数矩阵 (Hessian, 三阶导数, 参数偏导数)
    
    % 提取 Wilson-Cowan 参数
    c1 = par(1); c2 = par(2); c3 = par(3); c4 = par(4);
    P = par(5); Q = par(6);
    theta_u = par(7); theta_v = par(8);
    
    % 变量状态 (xx 里包含了时滞)
    u = xx(1,1);  % 当前时刻 u
    v = xx(2,1);  % 当前时刻 v
    u_tau1 = xx(1,2);  % 时滞 τ1 的 u
    v_tau2 = xx(2,3);  % 时滞 τ2 的 v
    
    % **Wilson-Cowan 的非线性激活函数**
    f_u = 1 ./ (1 + exp(-10 * (u - theta_u)));
    f_v = 1 ./ (1 + exp(-10 * (v - theta_v)));
    
    % **Wilson-Cowan 的一阶导数**
    dfdu = 10 * f_u .* (1 - f_u);
    dfdv = 10 * f_v .* (1 - f_v);
    
    % **Hessian (二阶导数)**
    if nx == 2 && np == 0
        d = zeros(2,2,2);  % Hessian 结构 (2×2×2)
        d(1,1,1) = -100 * f_u .* (1 - f_u) .* (1 - 2*f_u);
        d(2,2,2) = -100 * f_v .* (1 - f_v) .* (1 - 2*f_v);
        return;
    end
    
    % **三阶导数**
    if nx == 3 && np == 0
        d = zeros(2,2,2,2);  % 三阶导数结构 (2×2×2×2)
        d(1,1,1,1) = 1000 * f_u .* (1 - f_u) .* (1 - 6*f_u + 6*f_u.^2);
        d(2,2,2,2) = 1000 * f_v .* (1 - f_v) .* (1 - 6*f_v + 6*f_v.^2);
        return;
    end
    
    % **参数偏导数 (d²f/dx dP 和 d³f/dx² dP)**
    if np > 0
        d = zeros(2,2,np);
        if np == 1
            d(1,1,1) = dfdu; % 对 P 的导数
            d(2,2,1) = dfdv; % 对 Q 的导数
        end
        return;
    end
    
    % **如果 nx 不符合预期，则返回空**
    d = [];
end
