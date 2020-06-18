


%*********************************************************************
%*********************************************************************
%*********************************************************************
%功能： SICP算法的顶层实现，各部分程序实现见其M文件
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
xPre=310;           %激光扫描仪前一时刻位置横坐标
yPre=300;           %激光扫描仪前一时刻位置纵坐标
directionPre=0;     %激光扫描仪前一时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
xLate=330;          %激光扫描仪当前时刻位置横坐标
yLate=280;          %激光扫描仪当前时刻位置纵坐标
directionLate=20;    %激光扫描仪当前时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
numSample=1000;     %控制激光扫描仪采样点数，扫描一圈的采样点数为 samplingRange/180*samples
samplingRange=140;  %控制激光扫描仪采样范围，采样范围为(-samplingRange,samplingRange)
MapName='Map13.bmp';%选择实验地图（场景）
ExpData1FromCsv='Exp1_X_NoiseAdd_5dB.csv';%读取的csv文件名字/数据存为csv文件的名字
ExpData2FromCsv='Exp1_P_NoiseAdd_5dB.csv';%读取的csv文件名字/数据存为csv文件的名字
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

 %绘制数据点集平面图，可视化激光测距仪采集到的数据
figure(101);
title('添加高斯噪声后的激光数据点集');   %标题
xlabel('水平距离(cm)');                  %x轴
ylabel('竖直距离(cm)');                  %y轴
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off

%*********************************************************************
%*********************************************************************
%Step2：SICP迭代算法
%*********************************************************************
%Step2.1：相应参数设置
iteraNumber=100;                %设置算法迭代次数
p=0.4;
errorThreshold=0.5;             %设置提前终止迭代循环的误差阈值
Pi=P;                           %迭代第i次时的P点集，初始值设为原始P
Xi=X;
Po1=Pi;
setOut1=4:length(P(1,:));
distanceIndiv=zeros(size(P(1,:)));
FinalX=0;                   %记录横坐标位移
FinalY=0;                   %记录纵坐标位移
FinalTheta=0;               %记录旋转位移量
%*********************************************************************
%Step2.2：迭代循环
index=0;
for i0=1:iteraNumber
    %Step2.2.1：从X点集中寻找Pi点集中每个点的最近点,并得到对应点对
    setIn=[Pi;X];                       %两个2×n矩阵(Pi和X)构造一个4×n矩阵矩阵(setIn)，来满足ClosetPointMatch输入参数要求
    [setOut1]=ClosetPointMatch1(setIn);  %获取点集Pi的最近点集
    Xi=[setOut1(3,:);setOut1(4,:)];
    %Step2.2.2：计算X与Pi的误差 
    coordinateOffste=Xi-Pi;             %计算Xi与Pi的横纵坐标差值
    for j=1:length(coordinateOffste(1,:))
        distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%计算每个P点和X点之间的欧式距离
    end%for j=1:length(coordinateOffste(1,:)) 
    index=index+1;
    error(index)=sum(distanceIndiv)/length(distanceIndiv);%P集和X集的平均误差
    xMatrix(index)=FinalX;
    yMatrix(index)=FinalY;
    thetaMatrix(index)=FinalTheta;
    %Step2.2.3：判断是否提前终止循环，如果error<errorThreshold,则终止循环
    if error(index)<errorThreshold
        break;
    end
    
    for i11=1:40
    %Step2.2.4：根据得到的对应点对计算变换矩阵R、T和theta
    [R,T,theta]=CalculateTranMatrix1([Pi;Xi],p);
    FinalX=FinalX+T(1,1);                        %水平偏移量
    FinalY=FinalY+T(2,1);
    FinalTheta=FinalTheta+theta/pi*180;     %旋转偏移量
    Po1=Pi;%更新Po1,表示变换前的Pi
    for k=1:length(Pi(1,:))%更新Pi,表示变换后的Pi
        Pi(:,k)=R*Pi(:,k)+T;
    end
    %Step2.2.5：根据得到的R、T变换出P(i+1)
        diff1=Po1-Pi;
        for diff_i=1:size(diff1,2)%计算每列数据的內积的开方
            dotDiff(diff_i)=sqrt(dot(diff1(:,diff_i),diff1(:,diff_i)));
        end
        dotDiff=sort(dotDiff,'descend');%降序排序
        dual=dotDiff(1);%取最大值
        if(dual<1e-1)%小于误差阈值，则跳出当前循环
            break; 
        end
    end 
    
end%for i=1:iteraNumber
% FinalY=-FinalY;%定义的纵坐标水平向下，所以需要转换一下

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.3：求解MSE,任意一点的真实点为其欧氏距离最近点。

%Step2.3.1：从X点集中寻找Pi点集中每个点的最近点,并得到对应点对
PP=Pi;
XX=X;
setIn=[PP;XX];                       %两个2×n矩阵(Pi和X)构造一个4×n矩阵矩阵(setIn)，来满足ClosetPointMatch输入参数要求
[setOut1]=ClosetPointMatch1(setIn);  %获取点集Pi的最近点集
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.3.2：计算X与Pi的误差 
coordinateOffste=Xi-PP;             %计算Xi与Pi的横纵坐标差值
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%计算每个P点和X点之间的欧式距离
end%for j=1:length(coordinateOffste(1,:)) 
MSE=sum(distanceIndiv)/length(distanceIndiv);%P集和X集的平均误差

figure(103);
title('迭代误差');%标题
xlabel('迭代次数');%x轴
ylabel('误差迭代变化值');%y轴
hold on
for i=1:size(error,2)-1
    plot([i i+1],[error(i) error(i+1)],'r-x','MarkerSize',3,'LineWidth',3);  
end
hold off


figure(104);
title('坐标迭代变化值');%标题
xlabel('迭代次数');%x轴
ylabel('横/纵坐标迭代变化值');%y轴
hold on
for i=1:size(xMatrix,2)-1
    plot([i i+1],[xMatrix(i) xMatrix(i+1)],'g-x','MarkerSize',3,'LineWidth',3);  
    plot([i i+1],[yMatrix(i) yMatrix(i+1)],'m-x','MarkerSize',3,'LineWidth',3);
end
plot([1 size(xMatrix,2)],[xLate-xPre xLate-xPre],'b-','MarkerSize',3,'LineWidth',3);%横坐标收敛线
plot([1 size(yMatrix,2)],[yPre-yLate yPre-yLate],'r-','MarkerSize',3,'LineWidth',3);%纵坐标收敛线
hold off

figure(105)
title('旋转角度迭代变化值');%标题
xlabel('迭代次数');%x轴
ylabel('旋转角度迭代变化值');%y轴 
hold on
for i=1:size(thetaMatrix,2)-1
    plot([i i+1],[thetaMatrix(i) thetaMatrix(i+1)],'k-x','MarkerSize',3,'LineWidth',3);
end
plot([1 size(thetaMatrix,2)],[-directionLate -directionLate],'r-','LineWidth',3);%旋转角度收敛线

hold off



figure(106)
title('匹配效果图');%标题
xlabel('水平距离(cm)');%x轴
ylabel('竖直距离(cm)');%y轴
% set(gca,'xtick',[],'ytick',[]) % 同时去掉x、y轴的刻度
hold on
for i=1:size(Pi,2)
    plot([Pi(1,i) Pi(1,i)],[Pi(2,i) Pi(2,i)],'r.','MarkerSize',4);
end
for i=1:size(X,2)
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);
end
hold off

FinalTheta=-FinalTheta;




































