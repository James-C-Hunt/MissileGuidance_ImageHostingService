% ��2.1 ��άƽ��ս������-Ŀ�꽻սģ�ͷ���
clear
n=0;

% ״̬����ʼ��
VM = 3000.;     % �����ٶ�
VT = 1000.;     % Ŀ���ٶ�
XNT = 96.6;       % Ŀ����ٶ�
HEDEG = -20.;   % ָ�����ǣ��ȣ�
XNP = 4.;       % ��Ч��������ϵ��
RM1 = 0.;       % ������ʼλ��
RM2 = 10000.;
RT1 = 40000.;   % Ŀ���ʼλ��
RT2 = 10000.;
BETA=0.;        % Ŀ���ֵ�ٶȷ���
VT1=-VT*cos(BETA);  % Ŀ���ʼ�ٶ�
VT2=VT*sin(BETA);

HE=HEDEG/57.3;      % ���ǣ����ȣ�
T=0.;               % ʱ��
S=0.;               % ���ݼ�¼ʱ��

% ��Ŀ��ʼ���״̬������
RTM1=RT1-RM1;       %  ��Ŀ���� 
RTM2=RT2-RM2;       %  
RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
XLAM=atan2(RTM2,RTM1);              % ��Ŀ���߽�
XLEAD=asin(VT*sin(BETA+XLAM)/VM);   % ��ʼǰ�ý�
THET=XLAM+XLEAD;                    % ��ʼ�ٶ�ƫ��
VM1=VM*cos(THET+HE);                % ��ʼ�ٶ�
VM2=VM*sin(THET+HE);
VTM1 = VT1 - VM1;                   % ��Ŀ����ٶ�
VTM2 = VT2 - VM2;
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM;    % �ֽ��ٶ�

%% �����Ե�Ŀ��ϵ����
while VC >= 0       % ����ǰ
    % �趨���沽��
    if RTM < 1000   
        H=.0002;
    else
        H=.0002;
    end
    
    % ��¼״̬����ֵ
    BETAOLD=BETA;
    RT1OLD=RT1;
    RT2OLD=RT2;
    RM1OLD=RM1;
    RM2OLD=RM2;
    VM1OLD=VM1;
    VM2OLD=VM2;
    
    STEP=1;
    FLAG=0;
    while STEP <=1
        if FLAG==1
            STEP=2;
            %����״̬��Ԥ��
            BETA=BETA+H*BETAD;
            RT1=RT1+H*VT1;
            RT2=RT2+H*VT2;
            RM1=RM1+H*VM1;
            RM2=RM2+H*VM2;
            VM1=VM1+H*AM1;
            VM2=VM2+H*AM2;
            % ����ʱ�����
            T=T+H;
        end
        
        % ��Ŀ���ֵ����
        RTM1=RT1-RM1;           % ��Ŀ���λ��
        RTM2=RT2-RM2;
        RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
        VTM1=VT1-VM1;           % ��Ŀ����ٶ�
        VTM2=VT2-VM2;
        VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
        
        % ΢�ַ���
         % ������������ٶ�
        XLAM=atan2(RTM2,RTM1);  % ��Ŀ���߽�
        XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);  % ��Ŀ���߽��ٶ�
        XNC=XNP*VC*XLAMD;       % NVQ_ ���������õ����ٶ�ָ��
         % ������΢����
        AM1=-XNC*sin(XLAM);     % �������ٶ�
        AM2=XNC*cos(XLAM);
        VT1=-VT*cos(BETA);      % Ŀ���ٶ�
        VT2=VT*sin(BETA);
        BETAD=XNT/VT;           % Ŀ���ٶ�ƫ��
        FLAG=1;
    end
    FLAG=0;
    
    % ʹ�ö����������������״̬��
    BETA=.5*(BETAOLD+BETA+H*BETAD);
    RT1=.5*(RT1OLD+RT1+H*VT1);
    RT2=.5*(RT2OLD+RT2+H*VT2);
    RM1=.5*(RM1OLD+RM1+H*VM1);
    RM2=.5*(RM2OLD+RM2+H*VM2);
    VM1=.5*(VM1OLD+VM1+H*AM1);
    VM2=.5*(VM2OLD+VM2+H*AM2);
    
    % ״̬����¼
    S=S+H;
    if S >=.0009999
        S=0.;
        n=n+1;
        ArrayT(n)=T;
        ArrayRT1(n)=RT1;
        ArrayRT2(n)=RT2;
        ArrayRM1(n)=RM1;
        ArrayRM2(n)=RM2;
        ArrayXNCG(n)=XNC/32.2;
        ArrayRTM(n)=RTM;
    end
end
RTM

%% ͼ�����
figure(1)
plot(ArrayRT1,ArrayRT2,ArrayRM1,ArrayRM2),grid
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft) ')
ylabel('Altitude (Ft)')
hold on 
grid on

figure(2)
plot(ArrayT,ArrayXNCG),grid
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Time (sec)')
ylabel('Acceleration of missle (G)')
hold on
grid on

%% ���ݴ洢
clc
output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
save datfil.txt output -ascii
disp '*** Simulation Complete'
