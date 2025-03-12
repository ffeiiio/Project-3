#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 22:19:04 2025

@author: heyubo
"""

import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint  # 需要安装 `ddeint`，如果没有安装请先运行 `pip install ddeint`

# 定义 Wilson-Cowan 模型的时滞微分方程
def wc_dde(Y, t, d):
    # 读取当前状态
    u, v = Y(t)

    # 读取时滞状态（如果 t < tau，则使用初始函数）
    u_tau1, v_tau1 = Y(t - d['tau1']) if t > d['tau1'] else d['Y0']
    u_tau2, v_tau2 = Y(t - d['tau2']) if t > d['tau2'] else d['Y0']
    
    # 计算微分方程
    du_dt = -u + 1 / (1 + np.exp(-d['beta'] * (d['c1'] * u_tau1 + d['c2'] * v_tau2 + d['P'])))
    dv_dt = -v + 1 / (1 + np.exp(-d['beta'] * (d['c3'] * u_tau2 + d['c4'] * v_tau1 + d['Q'])))
    
    return np.array([du_dt, dv_dt])

# 定义参数
params = {
    'c1': 10, 'c2': -10, 'c3': 10, 'c4': 2,
    'P': -2.2, 'Q': -4, 'beta': 1,
    'tau1': 0.5, 'tau2': 0.5,  # 设定两个时滞
    'Y0': np.array([0.5, 3])  # 初始历史状态
}

# 定义初始历史函数
def init_history(t):
    return params['Y0']

# 计算时间演化
t_span = np.linspace(0, 50, 1000)
sol = ddeint(lambda Y, t: wc_dde(Y, t, params), init_history, t_span)

# 绘制结果
plt.figure(figsize=(8, 5))
plt.plot(t_span, sol[:, 0], label=r'$u(t)$', color='b')

plt.xlabel('Time')
plt.ylabel('Activity')
plt.title('Wilson-Cowan Model with Time Delays')
plt.legend()
ylim(0, 0.5)
plt.grid()
plt.show()