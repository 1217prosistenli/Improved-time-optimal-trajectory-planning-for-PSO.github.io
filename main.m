% 中北大学硕士论文采用3_5_3多段多项式插值求解
figure('name', '基于353的粒子群优化算法');
q0 = [1.0472 -0.9163 0.2618 -0.6545 -1.9722 0];   % 初始位置
q1 = [1.2357 -1.0420 0.1361 -0.8430 -1.7802 0];   % 中间节点1
q2 = [2.0595 -1.2305 -0.2409 -1.1572 -1.3352 0];   % 中间节点2
qf = [2.7053 -1.1676 -0.4922 -2.8972 -3.0718 0];  % 终止位置

N = 30; % 粒子群数量
kmax = 100;  % 迭代次数
D = 3;  % 维度为3，三段多项式

jj=1;
kk = 1;
retT = zeros(6,D);  % 记录所有关节的最终的最优位置

% 初始化种群
xx = 1.9*rand(N,D)+0.1;
vv = 2*rand(N,D)-1;
    

for m=1:6    % 每一个关节的求解
    % 用来计算传统PSO算法
    x = xx;
    v = vv;
    
    px = x;   % 个体最优位置
    pfit = ones(N,1)*inf;     % 粒子最优适应度
    gx = ones(1,D)*0.1;     % 全局最优位置
    gfit = inf;    % 全局最优适应度
    ret = zeros(1,kmax);  % 存储最优适应度
    
    % 用来计算改进后PSO算法
    x1 = xx;
    v1 = vv;
    px1 = x1;   % 个体最优位置
    pfit1 = ones(N,1)*inf;     % 粒子最优适应度
    gx1 = ones(1,D)*0.1;     % 全局最优位置
    gfit1 = inf;    % 全局最优适应度
    retx1 = zeros(N, D);  % 存储一个关节所有的时间 
    ret1 = zeros(1,kmax);  % 存储最优适应度
    
    k = 1;
    while k <= kmax
        [px, gx, gfit] = A(q0(m), q1(m), q2(m), qf(m), x, pfit, px, gx, gfit, N, D);   % 改进后的PSO算法计算个体最优位置以及全体最优位置
        [px1, gx1, gfit1] = A(q0(m), q1(m), q2(m), qf(m), x1, pfit1, px1, gx1, gfit1, N, D); % 传统的PSO算法计算个体最优位置以及全体最优位置
      
        [x1, v1] = i_pso(px1, gx1, x1, v1, k, kmax, N);  % 改进后的PSO算法更新位置与速度
        
        [x, v] = t_pso(px, gx, x, v, k, kmax, N);  % 传统的PSO算法更新位置与速度

        retx1(k,:) = gx1;    % 存储每一个关节最优位置变化
        ret1(k) = gfit1;       % 存储改进后的每一次迭代得到的最优适应度（也就是最佳的总时间）
        ret(k) = gfit;        % 存储传统的PSO每一次迭代得到的最优适应度（也就是最佳的总时间）
               
        k = k+1;    % 迭代次数增加
    end
    
    % 只绘制前三个关节的曲线变化
    if jj==m && m<4
        % 绘制每一个关节的各段时间曲线变化
        figure(kk);
        % retx1(:, i)表示的是第i段的所有时间
        plot(retx1(:, 1), 'r-');hold on;plot(retx1(:, 2), 'b--');hold on;plot(retx1(:, 3), 'm+');
        legend("第一段曲线时间","第二段曲线时间","第三段曲线时间");grid on;
        xlabel("迭代次数k");ylabel("t/s");title(["关节",m,"的时间变化"]);
        kk = kk+1;
        
        % 绘制每个关节的最优值gfit在迭代过程中的变化
        figure(kk);
        plot(ret1,'b--');hold on;
        plot(ret, 'r*');
        legend("改进PSO","PSO");xlabel("迭代次数k");ylabel("t/s");
        
        kk = kk+1;
        jj = jj+1;
    end
    
    retT(m, :) = gx;    % 存储所有关节的最优位置
end