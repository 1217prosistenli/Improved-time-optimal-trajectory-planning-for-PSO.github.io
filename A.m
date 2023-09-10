function [px, gx, gfit] = A(q0, q1, q2, qf, x, pfit, px, gx, gfit, N, D)  % 将最优位置和最优适应值返回
    Vmax = pi;  % 最大的速度值

    for i = 1:N     % 循环每一个粒子
        
        qd = [];    % 用来存储角速度

        % 计算系数a
        for j = 1:D  % 求解三个时间段的速度值
            if j==1
                % 第一段的三次多项式求解
                % 求解函数系数
                a10 = q0;
                a11 = 0;
                a12 = 0;
                a13 = (q1-q0)/x(i,j)^3;
                tz =0:0.01:x(i,j);
                % 求解角速度
                qz = a11+2*a12*tz+3*a13*tz.^2;
                qd = [qd, qz];
            elseif j==2
                % 第二段的五次多项式求解
                a20 = q1;
                a21 = 3*(q1-q0)/x(i,j-1);
                a22 = 3*(q1-q0)/x(i,j-1)^2;
                a23 = -(3*x(i,j)+12*x(i,j+1))*qf/(x(i,j)^2*x(i,j+1)^2)+(18*x(i,j-1)+9*x(i,j))*q0/(x(i,j-1)^2*x(i,j)^2)+(3*x(i,j)^2+12*x(i,j)*x(i,j+1)+10*x(i,j+1)^2)*q2/(x(i,j)^3*x(i,j+1)^2)-(10*x(i,j-1)^2+18*x(i,j-1)*x(i,j)+9*x(i,j)^2)*q1/(x(i,j-1)^2*x(i,j)^3);
                a24 = (6*x(i,j)+21*x(i,j+1))*qf/(x(i,j)^3*x(i,j+1)^2)-(24*x(i,j-1)+9*x(i,j))*q0/(x(i,j-1)^2*x(i,j)^3)-(6*x(i,j)^2+21*x(i,j)*x(i,j+1)+15*x(i,j+1)^2)*q2/(x(i,j)^4*x(i,j+1)^2)+(15*x(i,j-1)^2+24*x(i,j-1)*x(i,j)+9*x(i,j)^2)*q1/(x(i,j-1)^2*x(i,j)^4);
                a25 = (-3*x(i,j)-9*x(i,j+1))*qf/(x(i,j)^4*x(i,j+1)^2)+(9*x(i,j-1)+3*x(i,j))*q0/(x(i,j-1)^2*x(i,j)^4)+(3*x(i,j)^2+9*x(i,j)*x(i,j+1)+6*x(i,j+1)^2)*q2/(x(i,j)^5*x(i,j+1)^2)-(6*x(i,j-1)^2+9*x(i,j-1)*x(i,j)+3*x(i,j)^2)*q1/(x(i,j-1)^2*x(i,j)^5);
                tz = 0.01:0.01:x(i,j);
                % 求解角速度
                qz = a21+2*a22*tz+3*a23*tz.^2+4*a24*tz.^3+5*a25*tz.^4;
                qd = [qd, qz];
            else
                %第三段的三次多项式求解
                a30 = q2;
                a31 = 3*(qf-q2)/x(i,j);
                a32 = -3*(qf-q2)/x(i,j)^2;
                a33 = (qf-q2)/x(i,j)^3;
                tz = 0.01:0.01:x(i,j);
                % 求解角速度
                qz = a31+2*a32*tz+3*a33*tz.^2;
                qd = [qd, qz];
            end              
        end
        
        % 判断角速度是否超过最大值
        if max(abs(qd))>Vmax
            % 令其总时间为无穷大，便于后面进行去除
            pfit(i) = inf;
        else
            if sum(x(i,:)) < pfit(i)  % 比较个体适应度
                pfit(i) = sum(x(i,:)); % 更新个体历史最佳适应度
                px(i, :) = x(i, :); % 更新个体历史最佳位置
            end
        end

        % 找出所有粒子中最优的
        if pfit(i) < gfit     % 对比适应度
            gfit = pfit(i);     % 如果找到更优位置则进行更新最优值
            gx = px(i, :);  % 如果找到更优的位置则进行更新位置
        end
    end
end