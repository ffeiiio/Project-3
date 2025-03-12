#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:37:27 2025

@author: heyubo
"""

from pylab import *
from ddeint import ddeint

def wilson_cowan_model(Y,t,tau):
    u,v = Y(t)
    u_tau,v_tau = Y(t-tau)
    return array([
    -u + 1 / (1 + exp(-beta * (c1 * u_tau + c2 * v_tau + P))),
    -v + 1 / (1 + exp(-beta * (c3 * u_tau + c4 * v_tau + Q)))
])

# 参数
c1, c2, c3, c4 = 10, -10, 10, 2  # 连接强度
P, Q = -1.5, -4                        # 输入信号
beta = 1                               # Sigmoid斜率
tau=0

g = lambda t : array([0.5, 3])  # 初始历史状态
tt = linspace(0, 30, 10000)       # 求解时间范围


# 计算无时滞 (τ=0) 的解
yy = ddeint(wilson_cowan_model, g, tt, fargs=(tau,))
# 绘制相轨迹
plot(yy[:, 0], yy[:, 1], lw=2, label='No Delay (tau=0)')


xlim(0, 0.4)
ylim(0, 0.5)

legend() # display the legend
xlabel('u(t) - Excitatory Activity')
ylabel('v(t) - Inhibitory Activity')
title('Wilson-Cowan Phase Trajectory with Single Delay')
show()


