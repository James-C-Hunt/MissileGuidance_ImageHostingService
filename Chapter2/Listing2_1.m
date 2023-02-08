% 表2.1 二维平面战术导弹-目标交战模型仿真
clear
n=0;

% 状态量初始化
VM = 3000.;     % 导弹速度
VT = 1000.;     % 目标速度
XNT = 96.6;       % 目标加速度
HEDEG = -20.;   % 指向误差角（度）
XNP = 4.;       % 有效比例导引系数
RM1 = 0.;       % 导弹初始位置
RM2 = 10000.;
RT1 = 40000.;   % 目标初始位置
RT2 = 10000.;
BETA=0.;        % 目标初值速度方向
VT1=-VT*cos(BETA);  % 目标初始速度
VT2=VT*sin(BETA);

HE=HEDEG/57.3;      % 误差角（弧度）
T=0.;               % 时标
S=0.;               % 数据记录时标

% 弹目初始相对状态量计算
RTM1=RT1-RM1;       %  弹目距离 
RTM2=RT2-RM2;       %  
RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
XLAM=atan2(RTM2,RTM1);              % 弹目视线角
XLEAD=asin(VT*sin(BETA+XLAM)/VM);   % 初始前置角
THET=XLAM+XLEAD;                    % 初始速度偏角
VM1=VM*cos(THET+HE);                % 初始速度
VM2=VM*sin(THET+HE);
VTM1 = VT1 - VM1;                   % 弹目相对速度
VTM2 = VT2 - VM2;
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM;    % 抵近速度

%% 非线性弹目关系仿真
while VC >= 0       % 命中前
    % 设定仿真步长
    if RTM < 1000   
        H=.0002;
    else
        H=.0002;
    end
    
    % 记录状态量旧值
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
            %单步状态量预测
            BETA=BETA+H*BETAD;
            RT1=RT1+H*VT1;
            RT2=RT2+H*VT2;
            RM1=RM1+H*VM1;
            RM2=RM2+H*VM2;
            VM1=VM1+H*AM1;
            VM2=VM2+H*AM2;
            % 仿真时间更新
            T=T+H;
        end
        
        % 弹目相对值计算
        RTM1=RT1-RM1;           % 弹目相对位置
        RTM2=RT2-RM2;
        RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
        VTM1=VT1-VM1;           % 弹目相对速度
        VTM2=VT2-VM2;
        VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
        
        % 微分方程
         % 比例导引求加速度
        XLAM=atan2(RTM2,RTM1);  % 弹目视线角
        XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);  % 弹目视线角速度
        XNC=XNP*VC*XLAMD;       % NVQ_ 比例导引得到加速度指令
         % 求所需微分量
        AM1=-XNC*sin(XLAM);     % 导弹加速度
        AM2=XNC*cos(XLAM);
        VT1=-VT*cos(BETA);      % 目标速度
        VT2=VT*sin(BETA);
        BETAD=XNT/VT;           % 目标速度偏角
        FLAG=1;
    end
    FLAG=0;
    
    % 使用二阶龙格库塔法更新状态量
    BETA=.5*(BETAOLD+BETA+H*BETAD);
    RT1=.5*(RT1OLD+RT1+H*VT1);
    RT2=.5*(RT2OLD+RT2+H*VT2);
    RM1=.5*(RM1OLD+RM1+H*VM1);
    RM2=.5*(RM2OLD+RM2+H*VM2);
    VM1=.5*(VM1OLD+VM1+H*AM1);
    VM2=.5*(VM2OLD+VM2+H*AM2);
    
    % 状态量记录
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

%% 图像绘制
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

%% 数据存储
clc
output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
save datfil.txt output -ascii
disp '*** Simulation Complete'
