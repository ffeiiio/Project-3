clear;                           % clear variables
format compact
close all;                       % close figures
addpath('/Users/heyubo/Desktop/project 3/ddebiftool/DDE-Biftool-master/ddebiftool/');    % add ddebiftool folder to path
%#ok<*ASGLU,*NOPTS,*NASGU>

WC_tau = @() [8,9];  % Wilson-Cowan 的时滞参数 τ1=par(8), τ2=par(9)
ind_P = 5;           % P 在 par(5)，可以用来做 continuation
ind_taus = 8;    % Wilson-Cowan 的时滞参数索引


% 定义函数句柄
rhs_handle = @WC_sys_rhs;
deri_handle = @WC_sys_deri;

% 先用 set_funcs 创建 funcs 结构体
funcs = set_funcs(...
    'sys_rhs', rhs_handle, ...
    'sys_tau', @() [8,9], ...
    'sys_deri', deri_handle);

% 强制修正 sys_rhs
funcs = setfield(funcs, 'sys_rhs', rhs_handle);
% 强制修正 sys_deri
funcs = setfield(funcs, 'sys_deri', deri_handle);

% 提取参数（用于绘制分岔图）
getpar = @(x,i) arrayfun(@(p) p.parameter(i), x.point, 'UniformOutput', true);
% 提取状态变量（Wilson-Cowan 的 u 和 v）
getu = @(x) arrayfun(@(p) p.x(1), x.point, 'UniformOutput', true); % Wilson-Cowan 中 u 对应 x(1)
getv = @(x) arrayfun(@(p) p.x(2), x.point, 'UniformOutput', true); % Wilson-Cowan 中 v 对应 x(2)

% 提取 Hopf / Fold 分岔点
bgetpar = @(x,i,bif) arrayfun(@(p) p.parameter(i), x.point(br_getflags(x,bif)), 'UniformOutput', true);
bgetu = @(x,bif) arrayfun(@(p) p.x(1), x.point(br_getflags(x,bif)), 'UniformOutput', true);
bgetv = @(x,bif) arrayfun(@(p) p.x(2), x.point(br_getflags(x,bif)), 'UniformOutput', true);


%% 1️ 定义 Wilson-Cowan 模型的稳态解
% par: 参数向量 [c1, c2, c3, c4, P, Q, beta, tau1, tau2]
stst.kind = 'stst';  % 设定类型为 steady-state
stst.parameter = [10, -10, 10, 2, -2, -4, 1, 3.7, 1];  % 设定系统参数
stst.x = [0.5; 3];  % 初始稳态解


%% 2️ 线性稳定性分析 (Linear stability of initial equilibrium)
flag_newhheur = 1;  % 默认设置
method = df_mthod(funcs, 'stst', flag_newhheur);
method.stability.minimal_real_part = -1;  % 计算稳定性
[stst, success] = p_correc(funcs, stst, [], [], method.point);
stst.stability = p_stabil(funcs, stst, method.stability);

% 增加搜索更稳定特征值
method.stability.minimal_real_part = -10;
stst.stability = p_stabil(funcs, stst, method.stability);

%% 3️ 使用 SetupStst 自动构造稳态分支
contpar = ind_P;  % 继续参数 P
stst_branch = SetupStst(funcs, 'x', stst.x, ...
    'parameter', stst.parameter, ...
    'contpar', contpar, ...
    'max_step', [contpar, 0.1], ...
    'min_bound', [contpar, -6], ...
    'max_bound', [contpar,6]);

%% 4️ 继续计算分支
stst_branch0 = stst_branch;  % 创建 stst_branch0 作为初始分支

% 先计算正向方向的分支
[stst_branch0] = br_contn(funcs, stst_branch0, 100);

% 反转方向再计算
stst_branch0 = br_rvers(stst_branch0);
[stst_branch0] = br_contn(funcs, stst_branch0, 100);

%% 5️ 识别 Hopf / Fold 分岔点
[stst_branch_wbifs, stst_testfuncs] = LocateSpecialPoints(funcs, stst_branch0);
nunst_stst = GetStability(stst_branch_wbifs);


%分岔图u—P
%% 6️ 绘制 P vs u 的分岔图 (Bifurcation Diagram)
P_stst = getpar(stst_branch_wbifs, ind_P);  % 提取 P 参数
u_stst = getu(stst_branch_wbifs);  % Wilson-Cowan 模型中的 u

figure(2); clf;
ax2 = gca;
cla(ax2);

% 绘制稳定性区域（0 为稳定，>0 为不稳定）
plot(ax2, P_stst(nunst_stst == 0), u_stst(nunst_stst == 0), 'g.', ...  % 绿色表示稳定
     P_stst(nunst_stst > 0), u_stst(nunst_stst > 0), 'r.');  % 红色表示不稳定

hold on;
% 先检查 LocateSpecialPoints 是否找到 Hopf
hopf_indices = br_getflags(stst_branch_wbifs, 'hopf');
disp('Hopf Bifurcation Indices:');
disp(hopf_indices);

% 如果 hopf_indices 为空，手动计算 Hopf 位置
if isempty(hopf_indices)
    hopf_points = find((nunst_stst(1:end-1) == 0 & nunst_stst(2:end) > 0) | ...
                   (nunst_stst(1:end-1) > 0 & nunst_stst(2:end) == 0));


    disp('Detected Hopf Points:');
    disp(hopf_points);
else
    hopf_points = hopf_indices;
end

% 绘制 Hopf 分岔点
plot(P_stst(hopf_points), u_stst(hopf_points), 'mo', 'MarkerSize', 8, 'LineWidth', 2);
% 标注 Hopf / Fold 分岔点
plot(bgetpar(stst_branch_wbifs, ind_P, 'hopf'), bgetu(stst_branch_wbifs, 'hopf'), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
plot(bgetpar(stst_branch_wbifs, ind_P, 'fold'), bgetu(stst_branch_wbifs, 'fold'), 'mo', 'MarkerSize', 8, 'LineWidth', 2);

% 添加图例
stst_lgtext = {'stable', 'unstable', 'Hopf', 'Fold'};
legend(ax2, stst_lgtext, 'location', 'west');


% 轴标签 (Latex 格式)
xlabel('$P$', 'Interpreter', 'LaTeX');
ylabel('$u$', 'Interpreter', 'LaTeX');
