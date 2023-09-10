function [x, v] = i_pso(px, gx, x, v, k, kmax, N)
% �Ľ����PSO�㷨��λ�ú��ٶȵĸ���
        wmax = 0.9; % ����Ȩ��
        wmin = 0.4; % ��С��Ȩ��
        
     % ����Ȩ��
        if k <= 0.6*kmax
            w = (wmax-wmin)*(cos(k/kmax)).^2;
        else
            w = 0.2+0.1*rand();
        end
        
        % ѧϰ����c
        c1 = 2*(sin(k*pi/(2*kmax))).^2;
        c2 = 2- c1;
        
        % ����λ�ú��ٶ�
        v = w*v+c1*rand()*(px-x)+c2*rand()*(repmat(gx,N,1)-x);
        % Լ���ٶ�
        v(v<-1) = -1;
        v(v>1) = 1;
        
        x = x+v;        
        % Լ��ʱ��
        x(x<0.1) = 0.1;
        x(x>2) = 2;
end
