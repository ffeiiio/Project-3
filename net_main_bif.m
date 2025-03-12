clear;                           % 清除变量
format compact
close all;                       % 关闭所有图形窗口
addpath('/Users/heyubo/Desktop/project 3/ddebiftool/DDE-Biftool-master/ddebiftool/');    % 添加 DDE-Biftool 目录到路径
%#ok<*ASGLU,*NOPTS,*NASGU>

%% 1️⃣ 设置参数索引
ind_epsilon = 8;  % epsilon 在 par(8) 位置
WC_tau = @() [9,10];  % Wilson-Cowan 的时滞参数 τ=par(9), ρ=par(10)

%% 2️⃣ 定义函数句柄
rhs_handle = @net_sys_rhs;
deri_handle = @net_sys_deri;

funcs = set_funcs(...
    'sys_rhs', rhs_handle, ...
    'sys_tau', @() [9,10], ...  % 修正 τ 和 ρ 的索引
    'sys_deri', deri_handle);

% 强制修正 sys_rhs 和 sys_deri
funcs = setfield(funcs, 'sys_rhs', rhs_handle);
funcs = setfield(funcs, 'sys_deri', deri_handle);

% 提取参数
getpar = @(x,i) arrayfun(@(p) p.parameter(i), x.point, 'UniformOutput', true);
getu = @(x) arrayfun(@(p) p.x(1), x.point, 'UniformOutput', true); % 提取 u 变量

% 提取 Hopf / Fold 分岔点
bgetpar = @(x,i,bif) arrayfun(@(p) p.parameter(i), x.point(br_getflags(x,bif)), 'UniformOutput', true);
bgetu = @(x,bif) arrayfun(@(p) p.x(1), x.point(br_getflags(x,bif)), 'UniformOutput', true);

%% 3️⃣ 定义初始稳态解 (Steady-State)
% par: 参数向量 [c1, c2, c3, c4, P, Q, beta, epsilon, tau, rho]
stst.kind = 'stst';
stst.parameter = [-1, -0.4, -1, 0 , 0.5, 0.5, 60, 0.34, 0.5, 2.5];  
stst.x = [0.5; 0.5];  % 初始稳态解

%% 4️⃣ 线性稳定性分析
flag_newhheur = 1;
method = df_mthod(funcs, 'stst', flag_newhheur);
method.stability.minimal_real_part = -1;
[stst, success] = p_correc(funcs, stst, [], [], method.point);
stst.stability = p_stabil(funcs, stst, method.stability);


% 增加搜索更稳定特征值
method.stability.minimal_real_part = -10;
stst.stability = p_stabil(funcs, stst, method.stability);


%% 5️⃣ 使用 SetupStst 自动构造 `ϵ` 分支
contpar = ind_epsilon;  % 继续参数为 `ϵ`
stst_branch = SetupStst(funcs, 'x', stst.x, ...
    'parameter', stst.parameter, ...
    'contpar', contpar, ...
    'max_step', [contpar, 0.1], ...
    'min_bound', [contpar, -2], ...
    'max_bound', [contpar, 5]);

%% 6️⃣ 计算 `ϵ` 方向的分支
stst_branch0 = stst_branch;
[stst_branch0] = br_contn(funcs, stst_branch0, 100);  % 正向计算
stst_branch0 = br_rvers(stst_branch0);
[stst_branch0] = br_contn(funcs, stst_branch0, 100);  % 反向计算

%% 7️⃣ 识别 Hopf / Fold 分岔点
[stst_branch_wbifs, stst_testfuncs] = LocateSpecialPoints(funcs, stst_branch0);
nunst_stst = GetStability(stst_branch_wbifs);

%% 8️⃣ 绘制 `ϵ vs u` 分岔图
epsilon_stst = getpar(stst_branch_wbifs, ind_epsilon);  % 提取 ϵ 参数
u_stst = getu(stst_branch_wbifs);  % Wilson-Cowan 模型中的 u

figure(2); clf;
ax2 = gca;
cla(ax2);

% 绿色 = 稳定, 红色 = 不稳定
plot(ax2, epsilon_stst(nunst_stst == 0), u_stst(nunst_stst == 0), 'g.', ...
     epsilon_stst(nunst_stst > 0), u_stst(nunst_stst > 0), 'r.');

hold on;

% 先检查 LocateSpecialPoints 是否找到 Hopf
hopf_indices = br_getflags(stst_branch_wbifs, 'hopf');
disp('Hopf Bifurcation Indices:');
disp(hopf_indices);

% 如果 `hopf_indices` 为空，手动计算 Hopf 位置
if isempty(hopf_indices)
    hopf_points = find((nunst_stst(1:end-1) == 0 & nunst_stst(2:end) > 0) | ...
                   (nunst_stst(1:end-1) > 0 & nunst_stst(2:end) == 0));

    disp('Detected Hopf Points:');
    disp(hopf_points);
else
    hopf_points = hopf_indices;
end

% 绘制 Hopf 分岔点
plot(epsilon_stst(hopf_points), u_stst(hopf_points), 'mo', 'MarkerSize', 8, 'LineWidth', 2);

% 标注 Hopf / Fold 分岔点
plot(bgetpar(stst_branch_wbifs, ind_epsilon, 'hopf'), bgetu(stst_branch_wbifs, 'hopf'), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
plot(bgetpar(stst_branch_wbifs, ind_epsilon, 'fold'), bgetu(stst_branch_wbifs, 'fold'), 'mo', 'MarkerSize', 8, 'LineWidth', 2);

% 添加图例
stst_lgtext = {'stable', 'unstable', 'Hopf', 'Fold'};
legend(ax2, stst_lgtext, 'location', 'west');

% 轴标签 (Latex 格式)
xlabel('$\epsilon$', 'Interpreter', 'LaTeX');
ylabel('$u$', 'Interpreter', 'LaTeX');
title('Bifurcation Diagram: $\epsilon$ vs $u$', 'Interpreter', 'LaTeX');

