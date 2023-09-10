function [x, v] = t_pso(px, gx, x, v, k, kmax, N)
% 传统PSO算法对位置和速度的更新
    wmax = 0.9;
    wmin = 0.4;
    
    % 计算惯性权重
    w = wmax-(wmax-wmin)*(k/kmax);
    
    % 学习因子c
    c1 = 2;
    c2 = 2;
    
    % 更新位置和速度
    v = w*v+c1*rand()*(px-x)+c2*rand()*(repmat(gx,N,1)-x);
    % 约束速度
    v(v<-1) = -1;
    v(v>1) = 1;

    x = x+v;        
    % 约束时间
    x(x<0.1) = 0.1;
    x(x>2) = 2;
    
end