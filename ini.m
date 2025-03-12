clear; close all; clc;

% 添加 DDE-Biftool 目录到路径（确保路径正确）
addpath('/Users/heyubo/Desktop/project 3/ddebiftool/DDE-Biftool-master/ddebiftool/');

% 设置参数索引
ind_epsilon = 8;

% 定义函数句柄
rhs_handle = @net_sys_rhs;
deri_handle = @net_sys_deri;

funcs = set_funcs(...
    'sys_rhs', rhs_handle, ...
    'sys_tau', @() [9,10], ...  % 修正 τ 和 ρ 的索引
    'sys_deri', deri_handle);

% 设定 Wilson-Cowan 模型参数
% par: 参数向量 [c1, c2, c3, c4, P, Q, beta, epsilon, tau, rho]
base_par = [-1, -0.4, -1, 0 , 0.5, 0.5, 60, 0.34, 0.5, 2.5];

% 遍历范围
u_range = 0:0.01:1; % u 从 0 到 1，每次增加 0.05
v_range = 0:0.01:1; % v 从 0 到 1，每次增加 0.05

% 存储找到的稳态解
stable_solutions = [];

% 遍历 u 和 v
for u = u_range
    for v = v_range
        % 设定初值
        stst.kind = 'stst';
        stst.parameter = base_par;
        stst.x = [u; v];

        % 线性稳定性分析
        flag_newhheur = 1;
        method = df_mthod(funcs, 'stst', flag_newhheur);
        method.stability.minimal_real_part = -1;

        % 计算稳态解
        [stst, success] = p_correc(funcs, stst, [], [], method.point);
        
        % 如果成功找到平衡态解，则记录
        if success == 1
            stable_solutions = [stable_solutions; u, v]; %#ok<AGROW>
            fprintf('找到稳态解: u = %.2f, v = %.2f\n', u, v);
        end
    end
end

% 显示所有找到的稳态解
if isempty(stable_solutions)
    disp('未找到合适的初始稳态解');
else
    disp('找到的所有稳态解:');
    disp(stable_solutions);
end
