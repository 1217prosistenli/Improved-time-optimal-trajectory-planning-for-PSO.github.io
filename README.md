# Improved-time-optimal-trajectory-planning-for-PSO.github.io
基于时间最优的PSO轨迹规划
本文主要是基于353分段插值的轨迹规划，步骤如下：
1.通过PSO算法找到最优时间，PSO算法如下：
1.1 初始化粒子群，本文是寻找最优的时间，因此插值的每一段的时间t看成为粒子，t = 1.9*rand(N,D)+0.1,将初始化范围[0.1-2],其中N表示的是粒子群数量，D表示为维度，也就是三段时间，t1，t2，t3
for k = 1:kmax  % 迭代
  for i = 1:N  % 循环每一个粒子
    for j = 1:D  % 循环每个关节的每一段插值
      1.2 求解每一个粒子的角速度是否满足最大角速度要求，如果不满足的则将适应值设置为最大值（后期进行筛选时抛弃掉），满足的话就直接进行适应度求解，适应度函数f = min(t1+t2+t3)（在N个粒子中，三段时间之和最小那个）
        if sum(t1, t2, t3) < pbest(i)
          1.3 更新每一个个体最优值pbest(i) = sum(t1, t2, t3)
          1.4 并记录个体更优位置px(i) = (t1, t2, t3)
        if pbest(i) < gbest
          1.5 更新最优值gbest = pbest(i)
          1.6 并记录全局更优位置gx = px(i)
2.得到最优时间之后，采用353绘制出角位移、角速度以及角加速度的变化曲线
