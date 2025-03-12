function J = WC_sys_deri(xx, par, nx, np, v)
    J = [];

    % 🔹 提前定义 `f` 和 `f_prime`
    beta = par(7);
    f = @(x) 1 ./ (1 + exp(-beta * x));  % sigmoid
    f_prime = @(x) beta * f(x) .* (1 - f(x));  % sigmoid 的导数

    if length(nx)==1 && isempty(np) && isempty(v)
        % 计算 df/dx (对状态变量求偏导)
        f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
        f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));

        J = zeros(2,2);
        if any(nx == 0)  
            J(1,1) = -1;  
            J(2,2) = -1;
        end
        if any(nx == 1)  
            J(1,1) = par(1) * f1_x1;
            J(1,2) = par(2) * f1_x1;
        end
        if any(nx == 2)  
            J(2,1) = par(3) * f2_x2;
            J(2,2) = par(4) * f2_x2;
        end

    elseif isempty(nx) && length(np)==1 && isempty(v)  
        % 计算 df/dp (对参数求偏导)
        J = zeros(2,1);
        if np == 5 
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
            J(1,1) = f1_x1; 
            J(2,1) = 0; 
        elseif np == 8  
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
            f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));
            J(1,1) = -par(1) * f1_x1 * xx(1,2); 
            J(2,1) = -par(4) * f2_x2 * xx(2,2);
        elseif np == 9  
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
            f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));
            J(1,1) = -par(2) * f1_x1 * xx(2,3);  
            J(2,1) = -par(3) * f2_x2 * xx(1,3);
        end

    elseif length(nx)==1 && length(np)==1 && isempty(v)   
        J = zeros(2,1);
        if nx == 0  
            J(1,1) = -1;
            J(2,2) = -1;
        elseif nx == 1 && np == 5  
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));  
            J(1,1) = f1_x1; 
            J(2,1) = 0; 
        elseif nx == 2 && np == 5  
            f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));  
            J(1,1) = 0;  
            J(2,1) = f2_x2;  
        elseif nx == 1 && np == 8  
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
            f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));
            J(1,1) = -par(1) * f1_x1 * xx(1,2);  
            J(2,1) = -par(4) * f2_x2 * xx(2,2);  
        elseif nx == 2 && np == 9  
            f1_x1 = f_prime(par(1) * xx(1,2) + par(2) * xx(2,3) + par(5));
            f2_x2 = f_prime(par(3) * xx(1,3) + par(4) * xx(2,2) + par(6));
            J(1,1) = -par(2) * f1_x1 * xx(2,3);  
            J(2,1) = -par(3) * f2_x2 * xx(1,3);
        end

    elseif length(nx) == 2 && isempty(np) && ~isempty(v)  
        J = zeros(2);  
        if nx(1) == 1 && nx(2) == 1  
            f1_x1 = f_prime(xx(1,2));
            f2_x2 = f_prime(xx(2,2));
            J(1,1) = -2 * par(1) * f1_x1 * (1 - 2 * f1_x1) * v(1);
            J(2,2) = -2 * par(4) * f2_x2 * (1 - 2 * f2_x2) * v(2);
        elseif nx(1) == 2 && nx(2) == 2  
            f1_x1 = f_prime(xx(1,3));
            f2_x2 = f_prime(xx(2,3));
            J(1,1) = -2 * par(2) * f1_x1 * (1 - 2 * f1_x1) * v(1);
            J(2,2) = -2 * par(3) * f2_x2 * (1 - 2 * f2_x2) * v(2);
        elseif nx(1) == 1 && nx(2) == 2  
            f1_x1 = f_prime(xx(1,2));
            f2_x2 = f_prime(xx(2,3));
            J(1,2) = -par(2) * f1_x1 * (1 - 2 * f1_x1) * v(2);
            J(2,1) = -par(3) * f2_x2 * (1 - 2 * f2_x2) * v(1);
        end
    end

    disp('WC_sys_deri is running...');
    disp('par = ');
    disp(par);
end
