# Radar-Basic-Algorithm
一些可以在MATLAB中使用的基础雷达数据处理算法，目前可以分为以下三个部分：

## 回波数据处理
脉冲压缩、CA-CFAR恒虚警检测、PD测速、单脉冲测角等算法

## 滤波及数据融合
滤波包括了常增益滤波、卡尔曼(Kalman)滤波和扩展卡尔曼滤波(EKF)
数据融合采用BC和CC两种，基于KF和EKF实现

## 阵列信号处理
比较了不同设计准则下的输出信干噪比(OSINR)，包括MPDR、SMI、RMI、对角加载(DL)、子空间投影(SISP)
