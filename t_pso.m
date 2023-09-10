function [x, v] = t_pso(px, gx, x, v, k, kmax, N)
% ��ͳPSO�㷨��λ�ú��ٶȵĸ���
    wmax = 0.9;
    wmin = 0.4;
    
    % �������Ȩ��
    w = wmax-(wmax-wmin)*(k/kmax);
    
    % ѧϰ����c
    c1 = 2;
    c2 = 2;
    
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