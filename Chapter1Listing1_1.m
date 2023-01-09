close all
clc
T=0.;
S=0.;
Y=0.;% ״̬��1
YD=0.;% ״̬��2
X=1.;
H=.001;% ��ֵ���ֲ���
n=0.;
W=20;% rad/s
while T <=(1.-1e-5)
    YOLD=Y;             % ��¼��ǰ״̬��1
    YDOLD=YD;           % ��¼��ǰ״̬��2
    STEP=1;             % �ֲ����������־
    FLAG=0;             % �Ѿ�����״̬΢�ֱ�־
    while STEP <=1      % 
        if FLAG==1      % �Ѿ�����״̬΢��
            STEP=2;     % 
            Y=Y+H*YD;   % ����΢�ָ���״̬��1 
            YD=YD+H*YDD;% ����΢�ָ���״̬��2
            T=T+H;      % ����ʱ��
        end
        YDD=W*X-W*W*Y;	%΢�ַ��̣���״̬����״̬΢��
        FLAG=1;
    end
    FLAG=0;
    Y=.5*(YOLD+Y+H*YD);	   % ����Runge�CKutta���֣���ֵ+0.5tʱ�̹���+0.5t+hʱ�̹���
    YD=.5*(YDOLD+YD+H*YDD);% ����Runge�CKutta���֣���ֵ+0.5tʱ�̹���+0.5t+hʱ�̹���
    S=S+H;% ���ݴ������ڿ���
    if S >=.000999 % ���ݴ�������
        S=0.;
        n=n+1;
        ArrayT(n)=T;% ��¼��ǰʱ��
        ArrayY(n)=Y;% ��¼��ֵ�������ֵ
    end
end
figure

plot(ArrayT,ArrayY),grid  % ������ֵ���ֽ��
xlabel('Time (Sec)')
ylabel('y')

%% ���߼��룺�������ϱ任�Ƶ������۱�ʽ��
ArrayY_ = 1/W*(1-cos(W*ArrayT));
hold on
plot(ArrayT,ArrayY_,'--','LineWidth',2)
legend('Runge�CKutta','Closed-Form solution')
%%

clc
output=[ArrayT',ArrayY',ArrayY_'];   % ���߸��ģ�ͬʱ�����ʽ�����ֵ��
save datfil.txt output  -ascii
disp 'simulation finished'
