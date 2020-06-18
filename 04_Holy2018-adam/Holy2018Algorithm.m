%*********************************************************************
%*********************************************************************
%*********************************************************************
%功能： 主程序，各部分程序实现见其M文件
%作者：Shaofeng Wu 
%时间：2019.12.07
%邮箱：shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
clear
clc
%*********************************************************************
%*********************************************************************
%Step1：数据获取和数据预处理
%*********************************************************************
%Step1.1：设置激光扫描仪的前后两个位置，和采样点数
xPre=250;           %激光扫描仪前一时刻位置横坐标
yPre=250;           %激光扫描仪前一时刻位置纵坐标
directionPre=0;     %激光扫描仪前一时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
xLate=300;          %激光扫描仪当前时刻位置横坐标
yLate=300;          %激光扫描仪当前时刻位置纵坐标
directionLate=30;    %激光扫描仪当前时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
numSample=1000;     %控制激光扫描仪采样点数，扫描一圈的采样点数为 samplingRange/180*samples
samplingRange=140;  %控制激光扫描仪采样范围，采样范围为(-samplingRange,samplingRange)
MapName='Map13.bmp';%选择实验地图（场景）
ExpData1FromCsv='Exp14_X_NoiseAdd_5dB.csv';%读取的csv文件名字/数据存为csv文件的名字
ExpData2FromCsv='Exp14_P_NoiseAdd_5dB.csv';%读取的csv文件名字/数据存为csv文件的名字
DataSource=false;   %选择数据来源模式
%*********************************************************************
% Step1.2：获取激光扫描仪前后两个位置的坐标对应的X、P点集
if DataSource   %选择数据来源，true模式为重新生成含有新噪声的数据
                %选择数据来源，false模式为从对应的.csv文件中读取数据
                %因为添加的高斯噪声是随机的，每次生成的噪声不一样，
                %所以为了保持所有实验数据的一致性，需要将数据保存下来，
                %以方便做对比实验。
[setOutX0]=GetLaserDataPointSet(xPre,yPre,numSample,directionPre,11,MapName,samplingRange);   %初始位置设为（xPre,yPre），返回的数据点集记为模板点集X
[setOutP0]=GetLaserDataPointSet(xLate,yLate,numSample,directionLate,22,MapName,samplingRange);%移动位置设为（xLate,yLate），返回的数据点集记为待匹配点集P

X0=CoordinateTran(setOutX0);%坐标变换：极坐标变换为笛卡尔坐标
P0=CoordinateTran(setOutP0);%坐标变换：极坐标变换为笛卡尔坐标
 %可视化激光测距仪采集到的数据(加高斯噪声前)
figure(99);
title('原始激光数据点集');    %标题
hold on
for i=1:size(P0,2)
    plot([P0(1,i) P0(1,i)],[P0(2,i) P0(2,i)],'ro','MarkerSize',4); 
    plot([X0(1,i) X0(1,i)],[X0(2,i) X0(2,i)],'g.','MarkerSize',4);  
end
hold off
%*********************************************************************
% Step1.3：加入信噪比为5dB的高斯白噪声
[setOutX]=NoiseAddition(setOutX0,5);
[setOutP]=NoiseAddition(setOutP0,5);


%采样的数据存为.csv文件
%dlmwrite可以控制保存为CSV文件的数据精度位数
dlmwrite('Exp12_X_NoiseAdd_5dB11.csv', setOutX,'precision',32);
dlmwrite('Exp12_P_NoiseAdd_5dB11.csv', setOutP,'precision',32);
else
setOutX = csvread(ExpData1FromCsv);
setOutP = csvread(ExpData2FromCsv);
end %if true
%*********************************************************************
%Step1.4：坐标转换，将数据由极坐标系变换为笛卡尔坐标系
X=CoordinateTran(setOutX);%坐标变换：极坐标变换为笛卡尔坐标
P=CoordinateTran(setOutP);%坐标变换：极坐标变换为笛卡尔坐标 

%可视化激光测距仪采集到的数据(加了高斯噪声后)
figure(101);
title('添加高斯噪声后的激光数据点集');    %标题
% xlabel('水平距离(cm)');     %x轴
% ylabel('竖直距离(cm)');     %y轴
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'rx','MarkerSize',3); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'gx','MarkerSize',3);  
end
hold off

%*********************************************************************
%*********************************************************************
%Step2：Registration of lines算法
%*********************************************************************
%Step2.1：相应参数设置

%adam（自适应学习率算法）参数设置
parameters=[0;0;25];    %zeros(3,1);分别表示横坐标位移、纵坐标位移和旋转角度
learningRate = 0.1;     %初始学习率
b1 = 0.9;               %adam算法作者建议的默认值
b2 = 0.999;             %adam算法作者建议的默认值
e = 0.00000001;         %adam算法作者建议的默认值
mt=zeros(3,1);          %一阶动量
vt=zeros(3,1);          %二阶动量


thresholdError=100;     %阈值


Xk=X;
Pk=P;
% iter=0;
iter1=0;


for i=1:2500
%*********************************************************************
%Step2.2：调用函数计算直线梯度
[dx,dy,dTheta,sumError]=Registration_of_lines(Xk,Pk);
%*********************************************************************
%Step2.3：参数更新
%adam优化算法参数更新
gradient(1,1)=real(dx);       %横坐标位移梯度,只取实部
gradient(2,1)=real(dy);       %纵坐标位移梯度，只取实部
gradient(3,1)=real(dTheta);   %旋转角度梯度，只取实部

mt=b1*mt+(1-b1)*gradient;       %更新一阶动量
vt=b2*vt+(1-b2)*(gradient.^2);  %更新二阶动量
mtt=mt/(1-(b1.^(i+1)));         %误差修正
vtt=vt/(1-(b2.^(i+1)));         %误差修正
vtt_sqrt=[sqrt(vtt(1,1));sqrt(vtt(2,1));sqrt(vtt(3,1))];
%梯度下降法所求参数更新
parameters=parameters-learningRate*mtt./(vtt_sqrt+e);
finalDx=parameters(1,1);
finalDy=parameters(2,1);
finalTheta=parameters(3,1);
%*********************************************************************
%Step2.4：根据所得参数对点集进行相应变换
R=[cos(finalTheta/180*pi) -sin(finalTheta/180*pi);sin(finalTheta/180*pi) cos(finalTheta/180*pi)];       %旋转矩阵
T=[finalDx;finalDy];
for j=1:size(P,2)
        Pk(:,j)=R*P(:,j)+T;
end
%*********************************************************************
%Step2.5：记录误差、旋转角度和位移量在迭代中的变化（最后可视化迭代过程使用）
iter1=iter1+1;
res(1,iter1)=parameters(1,1);
res(2,iter1)=parameters(2,1);
res(3,iter1)=parameters(3,1);
error(iter1)=sumError; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.6:求解MSE,任意一点的真实点为其欧氏距离最近点。

%Step2.6.1：从X点集中寻找Pi点集中每个点的最近点,并得到对应点对
PP=Pk;
XX=Xk;
setIn=[PP;XX];                       %两个2×n矩阵(Pi和X)构造一个4×n矩阵矩阵(setIn)，来满足ClosetPointMatch输入参数要求
[setOut1]=ClosetPointMatch1(setIn);  %获取点集Pi的最近点集
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.6.2：计算X与Pi的误差 
coordinateOffste=Xi-PP;             %计算Xi与Pi的横纵坐标差值
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%计算每个P点和X点之间的欧式距离
end%for j=1:length(coordinateOffste(1,:)) 
MSE=sum(distanceIndiv)/length(distanceIndiv);%P集和X集的平均误差



%*********************************************************************
%*********************************************************************
%Step3：数据可视化
figure(13);
title('迭代误差');%标题
xlabel('迭代次数');%x轴
ylabel('误差迭代变化值');%y轴
hold on
for i=1:size(error,2)-1
    plot([i i+1],[error(i) error(i+1)],'k*-');
end
hold off

figure(14);
title('坐标迭代变化值');%标题
xlabel('迭代次数');%x轴
ylabel('横/纵坐标迭代变化值');%y轴
hold on
for i=1:size(res,2)-1
    plot([i i+1],[res(1,i) res(1,i+1)],'g*-');
    plot([i i+1],[xLate-xPre xLate-xPre],'b-','LineWidth',2);%收敛线
    
    plot([i i+1],[res(2,i) res(2,i+1)],'m*-');
    plot([i i+1],[yPre-yLate yPre-yLate],'r-','LineWidth',2);%收敛线
end
hold off

figure(15)
title('旋转角度迭代变化值');%标题
xlabel('迭代次数');%x轴
ylabel('旋转角度迭代变化值');%y轴 
hold on
for i=1:size(res,2)-1
    plot([i i+1],[res(3,i) res(3,i+1)],'k*-');
    plot([i i+1],[directionLate directionLate],'r-','LineWidth',2);%收敛线
end
hold off


figure(16)
title('匹配效果图');%标题
xlabel('水平距离(cm)');%x轴
ylabel('竖直距离(cm)');%y轴
hold on
for i=1:size(Pk,2)
    plot([Pk(1,i) Pk(1,i)],[Pk(2,i) Pk(2,i)],'r.','MarkerSize',4);
end
for i=1:size(Xk,2)
    plot([Xk(1,i) Xk(1,i)],[Xk(2,i) Xk(2,i)],'g.','MarkerSize',4);
end
hold off

































