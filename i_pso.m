function [x, v] = i_pso(px, gx, x, v, k, kmax, N)
% 改进后的PSO算法对位置和速度的更新
        wmax = 0.9; % 最大的权重
        wmin = 0.4; % 最小的权重
        
     % 计算权重
        if k <= 0.6*kmax
            w = (wmax-wmin)*(cos(k/kmax)).^2;
        else
            w = 0.2+0.1*rand();
        end
        
        % 学习因子c
        c1 = 2*(sin(k*pi/(2*kmax))).^2;
        c2 = 2- c1;
        
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
