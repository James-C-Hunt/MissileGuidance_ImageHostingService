close all
clc
T=0.;
S=0.;
Y=0.;% 状态量1
YD=0.;% 状态量2
X=1.;
H=.001;% 数值积分步长
n=0.;
W=20;% rad/s
while T <=(1.-1e-5)
    YOLD=Y;             % 记录当前状态量1
    YDOLD=YD;           % 记录当前状态量2
    STEP=1;             % 局部运算次数标志
    FLAG=0;             % 已经更新状态微分标志
    while STEP <=1      % 
        if FLAG==1      % 已经更新状态微分
            STEP=2;     % 
            Y=Y+H*YD;   % 利用微分更新状态量1 
            YD=YD+H*YDD;% 利用微分更新状态量2
            T=T+H;      % 更新时间
        end
        YDD=W*X-W*W*Y;	%微分方程：由状态量求状态微分
        FLAG=1;
    end
    FLAG=0;
    Y=.5*(YOLD+Y+H*YD);	   % 二阶Runge–Kutta积分：旧值+0.5t时刻估计+0.5t+h时刻估计
    YD=.5*(YDOLD+YD+H*YDD);% 二阶Runge–Kutta积分：旧值+0.5t时刻估计+0.5t+h时刻估计
    S=S+H;% 数据存数周期控制
    if S >=.000999 % 数据存数流程
        S=0.;
        n=n+1;
        ArrayT(n)=T;% 记录当前时刻
        ArrayY(n)=Y;% 记录数值积分输出值
    end
end
figure

plot(ArrayT,ArrayY),grid  % 绘制数值积分结果
xlabel('Time (Sec)')
ylabel('y')

%% 译者加入：绘制拉氏变换推导的理论闭式解
ArrayY_ = 1/W*(1-cos(W*ArrayT));
hold on
plot(ArrayT,ArrayY_,'--','LineWidth',2)
legend('Runge–Kutta','Closed-Form solution')
%%

clc
output=[ArrayT',ArrayY',ArrayY_'];   % 译者更改：同时保存闭式解和数值解
save datfil.txt output  -ascii
disp 'simulation finished'
