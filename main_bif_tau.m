clear;                           % clear variables
format compact
close all;                       % close figures
addpath('/Users/heyubo/Desktop/project 3/ddebiftool/DDE-Biftool-master/ddebiftool/');    % add ddebiftool folder to path
%#ok<*ASGLU,*NOPTS,*NASGU>

WC_tau = @() [8,9];  % Wilson-Cowan 的时滞参数 τ1=par(8), τ2=par(9)
ind_P = 5;           % P 在 par(5)，可以用来做 continuation
ind_taus = 8;        % Wilson-Cowan 的时滞参数索引（τ1）

% 定义函数句柄
rhs_handle = @WC_sys_rhs;
deri_handle = @WC_sys_deri;

% 先用 set_funcs 创建 funcs 结构体
funcs = set_funcs(...
    'sys_rhs', rhs_handle, ...
    'sys_tau', @() [8,9], ...
    'sys_deri', deri_handle);

% 让 τ2 (par(9)) 在整个计算过程中始终等于 τ1 (par(8))
funcs.get_comp = @(p) setfield(p, 'parameter', setfield(p.parameter, 9, p.parameter(8)));

% 提取参数（用于绘制分岔图）
getpar = @(x,i) arrayfun(@(p) p.parameter(i), x.point, 'UniformOutput', true);
% 提取状态变量（Wilson-Cowan 的 u 和 v）
getu = @(x) arrayfun(@(p) p.x(1), x.point, 'UniformOutput', true);
getv = @(x) arrayfun(@(p) p.x(2), x.point, 'UniformOutput', true);

% 提取 Hopf / Fold 分岔点
bgetpar = @(x,i,bif) arrayfun(@(p) p.parameter(i), x.point(br_getflags(x,bif)), 'UniformOutput', true);
bgetu = @(x,bif) arrayfun(@(p) p.x(1), x.point(br_getflags(x,bif)), 'UniformOutput', true);
bgetv = @(x,bif) arrayfun(@(p) p.x(2), x.point(br_getflags(x,bif)), 'UniformOutput', true);

%% 1️ 定义 Wilson-Cowan 模型的稳态解
stst.kind = 'stst';  % 设定类型为 steady-state
stst.parameter = [10, -10, 10, 2, -1.8, -4, 1, 3, 3];  % 设定系统参数
stst.x = [0.5; 3];  % 初始稳态解

%% 2️ 线性稳定性分析 (Linear stability of initial equilibrium)
flag_newhheur = 1;
method = df_mthod(funcs, 'stst', flag_newhheur);
method.stability.minimal_real_part = -1;
[stst, success] = p_correc(funcs, stst, [], [], method.point);
stst.stability = p_stabil(funcs, stst, method.stability);

% 增加搜索更稳定特征值
method.stability.minimal_real_part = -10;
stst.stability = p_stabil(funcs, stst, method.stability);

%% 7️ 重新构造以 τ1 = τ2 = τ 为 continuation 参数的稳态分支
contpar_tau = 8;  % 设定 τ1 (par(8)) 作为 continuation 参数

% 确保 τ2 也等于 τ1
stst.parameter(9) = stst.parameter(8);

stst_branch_tau = SetupStst(funcs, 'x', stst.x, ...
    'parameter', stst.parameter, ...
    'contpar', contpar_tau, ...
    'max_step', [contpar_tau, 0.001], ...
    'min_bound', [contpar_tau, 0], ...
    'max_bound', [contpar_tau, 3]);

%% 8️ 继续计算 τ 分支
stst_branch_tau = br_contn(funcs, stst_branch_tau, 300);
stst_branch_tau = br_rvers(stst_branch_tau);
stst_branch_tau = br_contn(funcs, stst_branch_tau, 300);

%% 9️ 识别 Hopf / Fold 分岔点
[stst_branch_tau_wbifs, stst_testfuncs_tau] = LocateSpecialPoints(funcs, stst_branch_tau);
nunst_stst_tau = GetStability(stst_branch_tau_wbifs);

%% 10️ 画 τ vs u 的分岔图
tau_stst = getpar(stst_branch_tau_wbifs, ind_taus);
u_stst = getu(stst_branch_tau_wbifs);

figure(3); clf;
ax3 = gca;
cla(ax3);

% 绘制稳定性区域（0 为稳定，>0 为不稳定）
plot(ax3, tau_stst(nunst_stst_tau == 0), u_stst(nunst_stst_tau == 0), 'g.', ...
     tau_stst(nunst_stst_tau > 0), u_stst(nunst_stst_tau > 0), 'r.');

hold on;

% 标注 Hopf / Fold 分岔点
if any(br_getflags(stst_branch_tau_wbifs, 'hopf'))
    plot(bgetpar(stst_branch_tau_wbifs, ind_taus, 'hopf'), ...
         bgetu(stst_branch_tau_wbifs, 'hopf'), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
end

if any(br_getflags(stst_branch_tau_wbifs, 'fold'))
    plot(bgetpar(stst_branch_tau_wbifs, ind_taus, 'fold'), ...
         bgetu(stst_branch_tau_wbifs, 'fold'), 'mo', 'MarkerSize', 8, 'LineWidth', 2);
end

% 添加图例
stst_lgtext_tau = {'stable', 'unstable', 'Hopf', 'Fold'};
legend(ax3, stst_lgtext_tau, 'location','south');

% 轴标签 (Latex 格式)
xlabel('$\tau$', 'Interpreter', 'LaTeX');
ylabel('$u$', 'Interpreter', 'LaTeX');
title('Bifurcation Diagram: $\tau$ vs $u$', 'Interpreter', 'LaTeX');
