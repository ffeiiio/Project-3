function J = net_sys_deri(xx, par, nx, np, v)
    J = [];

    % 🔹 提前定义 `f` 和 `f_prime`
    beta = par(7);
    f = @(x) 1 ./ (1 + exp(-beta * x));  % sigmoid
    f_prime = @(x) beta * f(x) .* (1 - f(x));  % sigmoid 的导数

    % 🔹 提取时滞变量
    u_tau = xx(1,2);  % u(t - tau)
    v_tau = xx(2,2);  % v(t - tau)
    u_rho = xx(1,3);  % u(t - rho)

    % 🔹 解析参数 (更新 rho 的索引)
    c1 = par(1); c2 = par(2); c3 = par(3); c4 = par(4);
    P = par(5); Q = par(6); 
    epsilon = par(8);  
    tau = par(9);      % ✅ 修正 τ
    rho = par(10);     % ✅ 修正 ρ（原本是 par(9)）

    % 🔹 计算 f'(x)
    f1_x1 = f_prime(c1 * u_tau + c2 * v_tau + epsilon * u_rho + P);  
    f2_x2 = f_prime(c3 * u_tau + c4 * v_tau + Q);

    % ========== 对状态变量求导 ==========
    if length(nx) == 1 && isempty(np) && isempty(v)
        J = zeros(2,2);
        if any(nx == 0)  % 对当前变量 x 求偏导
            J(1,1) = -1;  
            J(2,2) = -1;
        end
        if any(nx == 1)  % u 变量的偏导数
            J(1,1) = c1 * f1_x1;
            J(1,2) = c2 * f1_x1;
        end
        if any(nx == 2)  % v 变量的偏导数
            J(2,1) = c3 * f2_x2;
            J(2,2) = c4 * f2_x2;
        end

    % ========== 对参数 epsilon (ϵ) 求导 ==========
    elseif isempty(nx) && length(np) == 1 && isempty(v)  
        J = zeros(2,1);
        if np == 8  % ✅ ϵ 对 u 的贡献
            J(1,1) = f1_x1 * u_rho;  
            J(2,1) = 0;  
        end

        if np == 10  % ✅ ρ 对 u 的贡献 (原本是 np == 8，修正为 10)
            J(1,1) = -c3 * f1_x1 * u_rho;  
            J(2,1) = 0;
        end

    % ========== 对 (x, p) 混合求导 ==========
    elseif length(nx) == 1 && length(np) == 1 && isempty(v)   
        J = zeros(2,1);
        if nx == 1 && np == 8  % ✅ (du/dϵ)
            J(1,1) = f1_x1 * u_rho;  
            J(2,1) = 0;  
        end

        if nx == 1 && np == 10  % ✅ (du/dρ) (原本是 np == 8，修正为 10)
            J(1,1) = -c3 * f1_x1 * u_rho;  
            J(2,1) = 0;  
        end

    % ========== 二阶偏导 ==========
    elseif length(nx) == 2 && isempty(np) && ~isempty(v)  
        J = zeros(2,2);
        if nx(1) == 1 && nx(2) == 1  % d²u/dudϵ
            J(1,1) = c1 * f1_x1 * (1 - 2 * f1_x1) * v(1);
            J(2,2) = c3 * f2_x2 * (1 - 2 * f2_x2) * v(2);
        end
    end

    disp('✅ net_sys_deri 运行成功');
    disp('参数:');
    disp(par);
end
