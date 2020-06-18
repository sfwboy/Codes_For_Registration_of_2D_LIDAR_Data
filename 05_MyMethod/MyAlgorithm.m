%*********************************************************************
%*********************************************************************
%*********************************************************************
%功能：主程序，各部分程序实现见其M文件
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
directionLate=20;   %激光扫描仪当前时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
numSample=100;     %控制激光扫描仪采样点数，扫描一圈的采样点数为 samplingRange/180*samples
samplingRange=180;  %控制激光扫描仪采样范围，采样范围为(-samplingRange,samplingRange)
MapName='Map32.bmp';%选择实验地图（场景）
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

%可视化激光测距仪采集到的数据(加了高斯噪声后)
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
%Step2：激光测距仪数据的直线段提取
%*********************************************************************

%*********************************************************************
%Step2.1：数据集自适应分割
% [SetOut11]=SplitData(X,8,5);%数据点分割
% [SetOut22]=SplitData(P,8,5);%数据点分割
[SetOut11]=SplitData(X,2,5);%数据点分割
[SetOut22]=SplitData(P,2,5);%数据点分割
figure(1011);
hold on
for j=1:size(SetOut11,3)
    %plot(SetOut11(1,:,j),SetOut11(2,:,j),'o','Markersize',2,'Color',[1 0.9 0.5]);
    plot(SetOut11(1,:,j),SetOut11(2,:,j),'x','Markersize',3,'Color',[rand(1,1) rand(1,1) rand(1,1)]);
end
for j=1:size(SetOut22,3)
    %plot(SetOut11(1,:,j),SetOut11(2,:,j),'o','Markersize',2,'Color',[1 0.9 0.5]);
    plot(SetOut22(1,:,j),SetOut22(2,:,j),'o','Markersize',3,'Color',[rand(1,1) rand(1,1) rand(1,1)]);
end
    plot([0 0],[0 0],'r.','Markersize',20) %绘制原点（激光测距仪位置）
hold off


%*********************************************************************
%Step2.2：对分割后的数据子集，判断其是否为直线段特征，如果不是则再次分割

% [SetOut1]=SplitData1(SetOut11,15,3);%对非直线特征数据集只分割一次
% [SetOut2]=SplitData1(SetOut22,15,3);%对非直线特征数据集只分割一次
[SetOut1]=SplitData2(SetOut11,15,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
[SetOut2]=SplitData2(SetOut22,15,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集

%显示数据集分割和处理后的情况
figure(1012);
hold on
for j=1:size(SetOut1,3)
    %plot(SetOut11(1,:,j),SetOut11(2,:,j),'o','Markersize',2,'Color',[1 0.9 0.5]);
    plot(SetOut1(1,:,j),SetOut1(2,:,j),'x','Markersize',3,'Color',[rand(1,1) rand(1,1) rand(1,1)]);
end
for j=1:size(SetOut2,3)
    %plot(SetOut11(1,:,j),SetOut11(2,:,j),'o','Markersize',2,'Color',[1 0.9 0.5]);
    plot(SetOut2(1,:,j),SetOut2(2,:,j),'o','Markersize',3,'Color',[rand(1,1) rand(1,1) rand(1,1)]);
end
    plot([0 0],[0 0],'r.','Markersize',20)%绘制原点（激光测距仪位置）
hold off

%*********************************************************************
%Step2.3：直线段拟合
% [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1);%利用最小二乘法拟合直线，效果很好
% [Result_Ls2,NumFlag2]=LeastSquareLine1(SetOut2);%利用最小二乘法拟合直线，效果很好
% [Result_Ls1,NumFlag1]=LineFittingRANSAC(SetOut1);%利用RANSAC法拟合直线，效果一般
% [Result_Ls2,NumFlag2]=LineFittingRANSAC(SetOut2);%利用RANSAC法拟合直线，效果一般
[Result_Ls1,NumFlag1]=LineFitting(SetOut1);%RANSAC法改进为遍历数据集，拟合直线，效果很好
[Result_Ls2,NumFlag2]=LineFitting(SetOut2);%RANSAC法改进为遍历数据集，拟合直线，效果很好

%显示最终的拟合情况
LinePoint=zeros(2,2,size(SetOut1,3));%收集直线段的端点
figure(103);
hold on
for j=1:size(SetOut1,3)
    aa=rand(1,1);
    bb=rand(1,1);
    cc=rand(1,1);
    for i=1:size(SetOut1,2)
        if SetOut1(1,i,j)~=0 || SetOut1(2,i,j)~=0
%             plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'gx','MarkerSize',3); 
%             plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'x','Markersize',2,'Color',[aa bb cc]);
            plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'gx','MarkerSize',2);
        end
    end
end
for j=1:size(SetOut2,3)
    aa=rand(1,1);
    bb=rand(1,1);
    cc=rand(1,1);
    for i=1:size(SetOut2,2)
        if SetOut2(1,i,j)~=0 || SetOut2(2,i,j)~=0
            %plot([SetOut2(1,i,j) SetOut2(1,i,j)],[SetOut2(2,i,j) SetOut2(2,i,j)],'x','Markersize',2,'Color',[aa bb cc]);  
            plot([SetOut2(1,i,j) SetOut2(1,i,j)],[SetOut2(2,i,j) SetOut2(2,i,j)],'rx','MarkerSize',2);
        end
    end
end
%plot([0 0],[0 0],'r.','Markersize',20);
for j=1:size(SetOut1,3)
    if abs(SetOut1(1,1,j)-SetOut1(1,NumFlag1(j),j))>5
        plot([SetOut1(1,1,j) SetOut1(1,NumFlag1(j),j)],...
        [Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j) ...
        Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j)],'b-v','LineWidth',1,'MarkerSize',3);
        LinePoint(1,1,j)=SetOut1(1,1,j);%起始点横坐标
        LinePoint(1,2,j)=SetOut1(1,NumFlag1(j),j);%结束点横坐标
        LinePoint(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%起始点纵坐标
        LinePoint(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%结束点纵坐标
    else
        plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
        (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
        [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-v','LineWidth',1,'MarkerSize',3);
        LinePoint(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%起始点横坐标
        LinePoint(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%结束点横坐标
        LinePoint(2,1,j)=SetOut1(2,1,j);%起始点纵坐标
        LinePoint(2,2,j)=SetOut1(2,NumFlag1(j),j);%结束点纵坐标
    end 
end
for j=1:size(SetOut2,3)
    if abs(SetOut2(1,1,j)-SetOut2(1,NumFlag2(j),j))>5
        plot([SetOut2(1,1,j) SetOut2(1,NumFlag2(j),j)],...
        [Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j) ...
        Result_Ls2(1,j)*SetOut2(1,NumFlag2(j),j)+Result_Ls2(2,j)],'k-o','LineWidth',1,'MarkerSize',3);
    else
        plot([(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j) ...
        (SetOut2(2,NumFlag2(j),j)-Result_Ls2(2,j))/Result_Ls2(1,j)],...
        [SetOut2(2,1,j) SetOut2(2,NumFlag2(j),j)],'k-o','LineWidth',1,'MarkerSize',3);
    end 
end
hold off





%*********************************************************************
%*********************************************************************
%Step3：计算直线段旋转角度
%*********************************************************************

%*********************************************************************
%Step3.1：根据直线参数y=ax+b，计算直线倾斜角度
for i=1:size(Result_Ls1,2)
    dist1(i)=abs(Result_Ls1(2,i)/sqrt(Result_Ls1(1,i)^2+1));
    theta1(i)=atan(Result_Ls1(1,i));
end
for i=1:size(Result_Ls2,2)
    dist2(i)=abs(Result_Ls2(2,i)/sqrt(Result_Ls2(1,i)^2+1));
    theta2(i)=atan(Result_Ls2(1,i));
end    
%*********************************************************************
%Step3.2：角度作差，顺便将角度由弧度制转为角度制
DiffTheta=0;
AddTheta=1;
for i=1:size(theta1,2)
    for j=1:size(theta2,2)
        DiffTheta(1,AddTheta)=theta1(i)/pi*180-theta2(j)/pi*180;
%         DiffTheta(2,AddTheta)=theta1(i)/pi*180-theta2(j)/pi*180;
        DiffTheta(2,AddTheta)=DiffTheta(1,AddTheta);
        AddTheta=AddTheta+1;
    end
end
%*********************************************************************
%Step3.3：使用密度聚类来计算角度差
figure(105);
hold on
[DbscanNum,DbscanCell] = DBSCAN(DiffTheta,1,2);%密度聚类
hold off
%选取元素最多的簇为关键簇
SizeCluster=zeros(1,2);
tempCell=DbscanCell{1,1};
SizeCluster(1,1)=size(tempCell,2);%记录当前簇包含的数据个数
SizeCluster(1,2)=1;               %记录对应的标记
%搜索关键簇
for i=2:DbscanNum
    tempCell=DbscanCell{i,1};
    if size(tempCell,2)>SizeCluster(1,1)%如果当前簇元素大于前一个，则更新
        SizeCluster(1,1)=size(tempCell,2);%记录当前簇包含的数据个数
        SizeCluster(1,2)=i;%记录对应的标记
    end
end
%导出关键簇
FlagDBC=DbscanCell{SizeCluster(1,2),1};
FinalTheta=0;
%对关键簇的所有样本加权平均即可得到最后的旋转角度
for i=1:size(FlagDBC,2)
    FinalTheta=FinalTheta+DiffTheta(1,FlagDBC(1,i));
end
FinalTheta=FinalTheta/size(FlagDBC,2);
% FinalTheta(2,1)=FinalTheta(2,1)/size(FlagDBC,2);

%*********************************************************************
%Step3.4：将P点集所有点根据得到的角度旋转回来,为下面的计算位移做准备
%得到旋转角度
TranMatrix=[cos(-FinalTheta/180*pi) -sin(-FinalTheta/180*pi);...
            sin(-FinalTheta/180*pi) cos(-FinalTheta/180*pi)];
SetOut2Tra=SetOut2;
%将P点集所有点旋转
for j=1:size(SetOut2,3)
    PP=SetOut2(:,:,j)';
    for i=1:size(PP,1)
        P1=PP(i,:)*TranMatrix;
        SetOut2Tra(1,i,j)=P1(1,1);
        SetOut2Tra(2,i,j)=P1(1,2);
    end
end
figure(104);
hold on
%绘制出旋转之后的直线
for j=1:size(SetOut1,3)
    for i=1:size(SetOut1,2)
        if SetOut1(1,i,j)~=0 || SetOut1(2,i,j)~=0
            plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'gx','MarkerSize',2);    
        end
    end
end
for j=1:size(SetOut2Tra,3)
    for i=1:size(SetOut2Tra,2)
        if SetOut2Tra(1,i,j)~=0 || SetOut2Tra(2,i,j)~=0
            plot([SetOut2Tra(1,i,j) SetOut2Tra(1,i,j)],[SetOut2Tra(2,i,j) SetOut2Tra(2,i,j)],'rx','MarkerSize',2);    
        end
    end
end
hold off

%利用数学关系根据旋转前的直线参数计算出旋转后的直线参数，其实也可以对旋转后的点集
%拟合直线，但是为了保持直线参数的一致性，所以采用数学关系转换参数
Line2StartPoint=zeros(2,size(SetOut2,3));
for j=1:size(SetOut2,3)
    if abs(SetOut2(1,1,j)-SetOut2(1,NumFlag2(j),j))>5
        %每条直线选取起始点，用来进行下面的直线方程求解
        Line2StartPoint(1,j)=SetOut2(1,1,j);%起始点横坐标
        Line2StartPoint(2,j)=Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j);%起始点纵坐标
    else
        %每条直线选取起始点，用来进行下面的直线方程求解
        Line2StartPoint(1,j)=(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j);%起始点横坐标
        Line2StartPoint(2,j)=SetOut2(2,1,j);%起始点纵坐标
    end 
end

%已经知道旋转后直线的斜率Result_Ls2Tra和旋转矩阵TranMatrix1
%以及直线上的一个起始点Line2StartPoint
%所以可以求解旋转前直线的斜率K和截距
Result_Ls2Tra=zeros(2,size(Result_Ls2,2));
for i=1:size(Result_Ls2Tra,2)
     Result_Ls2Tra(1,i)=tan(theta2(i)+FinalTheta/180*pi);
end
for i=1:size(Result_Ls2,2)
    temp=Line2StartPoint(:,i)'*(TranMatrix);%将直线上的点旋转
    Result_Ls2Tra(2,i)=temp(2)-Result_Ls2Tra(1,i)*temp(1);   
end
figure(104);
hold on
%显示旋转后的直线拟合情况
Line1Point=zeros(2,2,size(SetOut1,3));%收集直线段的端点
Line2TraPoint=zeros(2,2,size(SetOut2Tra,3));%收集直线段的端点
for j=1:size(SetOut1,3)
    if abs(SetOut1(1,1,j)-SetOut1(1,NumFlag1(j),j))>5
        plot([SetOut1(1,1,j) SetOut1(1,NumFlag1(j),j)],...
        [Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j) ...
        Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j)],'b-v','LineWidth',1,'MarkerSize',3);
        Line1Point(1,1,j)=SetOut1(1,1,j);%起始点横坐标
        Line1Point(1,2,j)=SetOut1(1,NumFlag1(j),j);%结束点横坐标
        Line1Point(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%起始点纵坐标
        Line1Point(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%结束点纵坐标
    else
        plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
        (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
        [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-v','LineWidth',1,'MarkerSize',3);
        Line1Point(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%起始点横坐标
        Line1Point(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%结束点横坐标
        Line1Point(2,1,j)=SetOut1(2,1,j);%起始点纵坐标
        Line1Point(2,2,j)=SetOut1(2,NumFlag1(j),j);%结束点纵坐标
    end 
end
for j=1:size(SetOut2Tra,3)
    if abs(SetOut2Tra(1,1,j)-SetOut2Tra(1,NumFlag2(j),j))>5
        plot([SetOut2Tra(1,1,j) SetOut2Tra(1,NumFlag2(j),j)],...
        [Result_Ls2Tra(1,j)*SetOut2Tra(1,1,j)+Result_Ls2Tra(2,j) ...
        Result_Ls2Tra(1,j)*SetOut2Tra(1,NumFlag2(j),j)+Result_Ls2Tra(2,j)],'k-o','LineWidth',1,'MarkerSize',3);
        Line2TraPoint(1,1,j)=SetOut2Tra(1,1,j);%起始点横坐标
        Line2TraPoint(1,2,j)=SetOut2Tra(1,NumFlag2(j),j);%结束点横坐标
        Line2TraPoint(2,1,j)=Result_Ls2Tra(1,j)*SetOut2Tra(1,1,j)+Result_Ls2Tra(2,j);%起始点纵坐标
        Line2TraPoint(2,2,j)=Result_Ls2Tra(1,j)*SetOut2Tra(1,NumFlag2(j),j)+Result_Ls2Tra(2,j);%结束点纵坐标
    else
        plot([(SetOut2Tra(2,1,j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j) ...
        (SetOut2Tra(2,NumFlag2(j),j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j)],...
        [SetOut2Tra(2,1,j) SetOut2Tra(2,NumFlag2(j),j)],'k-o','LineWidth',1,'MarkerSize',3);
        Line2TraPoint(1,1,j)=(SetOut2Tra(2,1,j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j);%起始点横坐标
        Line2TraPoint(1,2,j)=(SetOut2Tra(2,NumFlag2(j),j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j);%结束点横坐标
        Line2TraPoint(2,1,j)=SetOut2Tra(2,1,j);%起始点纵坐标
        Line2TraPoint(2,2,j)=SetOut2Tra(2,NumFlag2(j),j);%结束点纵坐标
    end 
end
hold off
for i=1:size(Result_Ls2Tra,2)
%     dist2Tra(i)=abs(Result_Ls2Tra(2,i)/sqrt(Result_Ls2Tra(1,i)^2+1));%原点到直线距离
    dist2Tra(i)=dist2(i);%因为直线旋转只改变θ，不改变ρ，所以dist2Tra(i)=abs(Result_Ls2Tra(2,i)/sqrt(Result_Ls2Tra(1,i)^2+1))也正确，两者理论上等效（实际近乎相等）
    theta2Tra(i)=atan(Result_Ls2Tra(1,i));
end 


%*********************************************************************
%*********************************************************************
%Step4：计算直线段位移
%*********************************************************************


%*********************************************************************
%Step4.1：计算直线间的距离，因为理论上大部分直线都是不平行的（最多近似平行）
%         所以直线距离计算原理是：取一条直线段的中点，然后计算该中点到另外直
%         线段的距离。

%计算每条线段的中点
Line1CenterPoints=zeros(2,size(Line1Point,3));
Line2TraCenterPoints=zeros(2,size(Line2TraPoint,3));
for i=1:size(Line1CenterPoints,2)
        startX=Line1Point(1,1,i);
        startY=Line1Point(2,1,i);
        endX=Line1Point(1,2,i);
        endY=Line1Point(2,2,i);
        xMid=startX-(startX-endX)/2;
        yMid=startY-(startY-endY)/2;
        Line1CenterPoints(1,i)=xMid;
        Line1CenterPoints(2,i)=yMid;
end
for i=1:size(Line2TraCenterPoints,2)
        startX=Line2TraPoint(1,1,i);
        startY=Line2TraPoint(2,1,i);
        endX=Line2TraPoint(1,2,i);
        endY=Line2TraPoint(2,2,i);
        xMid=startX-(startX-endX)/2;
        yMid=startY-(startY-endY)/2;
        Line2TraCenterPoints(1,i)=xMid;
        Line2TraCenterPoints(2,i)=yMid;
end

%候选直线段对
Add=1;
CorrectLinePairsFlag=0;
for i=1:size(theta1,2) 
    for j=1:size(theta2Tra,2)
        if abs(theta1(i)-theta2Tra(j))<2/180*pi%2度
            CorrectLinePairsFlag(1,Add)=i;
            CorrectLinePairsFlag(2,Add)=j;
            Add=Add+1;
        end
    end
end




%根据直线对计算距离
DistLines=zeros(1,size(CorrectLinePairsFlag,2));
for i=1:size(CorrectLinePairsFlag,2)
    xMid=Line1CenterPoints(1,CorrectLinePairsFlag(1,i));
    yMid=Line1CenterPoints(2,CorrectLinePairsFlag(1,i));
    DistLines(i)=abs((Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))*xMid-yMid+Result_Ls2Tra(2,CorrectLinePairsFlag(2,i)))/sqrt(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
end


%*********************************************************************
%Step4.2：角度校正后，点集P的直线段的法向量


%角度校正后，点集P的直线段的法向量
VectorN=zeros(2,size(CorrectLinePairsFlag,2));
for i=1:size(CorrectLinePairsFlag,2)
    if abs(Result_Ls1(1,CorrectLinePairsFlag(1,i)))> 12 || abs(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i)))>12
        x1=-Result_Ls1(2,CorrectLinePairsFlag(1,i))/Result_Ls1(1,CorrectLinePairsFlag(1,i));
        x2=-Result_Ls2Tra(2,CorrectLinePairsFlag(2,i))/Result_Ls2Tra(1,CorrectLinePairsFlag(2,i));
        a=1;
        if x1>x2
            VectorN(1,i)=-sqrt(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2/(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
            VectorN(2,i)=-VectorN(1,i)/Result_Ls2Tra(1,CorrectLinePairsFlag(2,i));
        else
            VectorN(1,i)=sqrt(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2/(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
            VectorN(2,i)=-VectorN(1,i)/Result_Ls2Tra(1,CorrectLinePairsFlag(2,i));
        end
    else
        if Result_Ls1(2,CorrectLinePairsFlag(1,i))>Result_Ls2Tra(2,CorrectLinePairsFlag(2,i))
            VectorN(1,i)=Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))*sqrt(1/(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
            VectorN(2,i)=-VectorN(1,i)/Result_Ls2Tra(1,CorrectLinePairsFlag(2,i));
        else
            VectorN(1,i)=-Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))*sqrt(1/(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
            VectorN(2,i)=-VectorN(1,i)/Result_Ls2Tra(1,CorrectLinePairsFlag(2,i));
        end
    end
end

%*********************************************************************
%Step4.3：将所有的直线对分割成两对直线段为一组的大集合
PreValue=CorrectLinePairsFlag(1,1);
LineMatrix(1,1,1)=CorrectLinePairsFlag(1,1);
LineMatrix(2,1,1)=CorrectLinePairsFlag(2,1);
LineMatrix(3,1,1)=-VectorN(1,1);
LineMatrix(4,1,1)=-VectorN(2,1);
LineMatrix(5,1,1)=DistLines(1);
add11=1;
add12=1;
for i=2:size(CorrectLinePairsFlag,2)
    if CorrectLinePairsFlag(1,i)~=PreValue
        PreValue=CorrectLinePairsFlag(1,i);
        add11=add11+1;
        add12=1;
        LineMatrix(1,1,add11)=CorrectLinePairsFlag(1,i);
        LineMatrix(2,1,add11)=CorrectLinePairsFlag(2,i);
        LineMatrix(3,1,add11)=-VectorN(1,i);
        LineMatrix(4,1,add11)=-VectorN(2,i);
        LineMatrix(5,1,add11)=DistLines(i);
    else
        add12=1+add12;
        LineMatrix(1,add12,add11)=CorrectLinePairsFlag(1,i);
        LineMatrix(2,add12,add11)=CorrectLinePairsFlag(2,i);
        LineMatrix(3,add12,add11)=-VectorN(1,i);
        LineMatrix(4,add12,add11)=-VectorN(2,i);  
        LineMatrix(5,add12,add11)=DistLines(i);
    end
end
addFinal=1;
for i=1:size(LineMatrix,3)-1
    for j=i+1:size(LineMatrix,3)
        if abs(theta1(LineMatrix(1,1,i))-theta1(LineMatrix(1,1,j)))<5/180*pi%两条直线的夹角大于设定的5度，避免得到的两对线段对为平行直线
            continue;%如果夹角小于阈值，跳过当前集合
        end 
        temp1=LineMatrix(:,:,i);
        temp2=LineMatrix(:,:,j);
        for m=1:size(temp1,2)
            if temp1(1,m)==0
                break;
            end
            for n=1:size(temp2,2)
                if temp2(1,n)==0
                    break;
                end
                FinalLinePair(1,addFinal)=temp1(1,m);
                FinalLinePair(2,addFinal)=temp1(2,m);
                FinalLinePair(3,addFinal)=temp1(3,m);
                FinalLinePair(4,addFinal)=temp1(4,m);
                FinalLinePair(5,addFinal)=temp1(5,m);
            
                FinalLinePair(6,addFinal)=temp2(1,n);
                FinalLinePair(7,addFinal)=temp2(2,n);
                FinalLinePair(8,addFinal)=temp2(3,n);
                FinalLinePair(9,addFinal)=temp2(4,n);
                FinalLinePair(10,addFinal)=temp2(5,n);
                addFinal=addFinal+1;
            end
        end
    end
end

%*********************************************************************
%Step4.4：最小二乘法拟合位移

for i=1:size(FinalLinePair,2)
    X_Ls(1,1)=FinalLinePair(3,i);
    X_Ls(1,2)=FinalLinePair(4,i);
    Y_Ls(1,1)=FinalLinePair(5,i);
    
    X_Ls(2,1)=FinalLinePair(8,i);
    X_Ls(2,2)=FinalLinePair(9,i);
    Y_Ls(2,1)=FinalLinePair(10,i);
    
    FinalVectorMatrix(:,i)=inv(X_Ls'*X_Ls)*X_Ls'*Y_Ls;
end


%*********************************************************************
%Step4.5：使用密度聚类来优化位移
figure(106);
hold on
[DbscanNum1,DbscanCell1] = DBSCAN(FinalVectorMatrix,1.5,2);
hold off
%选取元素最多的簇为关键簇
SizeCluster1=zeros(1,2);
tempCell1=DbscanCell1{1,1};
SizeCluster1(1,1)=size(tempCell1,2);%记录当前簇包含的数据个数
SizeCluster1(1,2)=1;                %记录对应的标记
%搜索关键簇
for i=2:DbscanNum1
    tempCell1=DbscanCell1{i,1};
    if size(tempCell1,2)>SizeCluster1(1,1)%如果当前簇元素大于前一个，则更新
        SizeCluster1(1,1)=size(tempCell1,2);%记录当前簇包含的数据个数
        SizeCluster1(1,2)=i;%记录对应的标记
    end
end
%导出关键簇
FlagDBC1=DbscanCell1{SizeCluster1(1,2),1};
FinalVector=zeros(2,1);
%对关键簇的所有样本加权平均即可得到最后的旋转角度
for i=1:size(FlagDBC1,2)
    FinalVector(1,1)=FinalVector(1,1)+FinalVectorMatrix(1,FlagDBC1(1,i));
    FinalVector(2,1)=FinalVector(2,1)+FinalVectorMatrix(2,FlagDBC1(1,i));
end
FinalVector(1,1)=FinalVector(1,1)/size(FlagDBC1,2);
FinalVector(2,1)=FinalVector(2,1)/size(FlagDBC1,2);


figure(107)
title('匹配效果图');%标题
xlabel('水平距离(cm)');%x轴
ylabel('竖直距离(cm)');%y轴
%set(gca,'xtick',[],'ytick',[]) % 同时去掉x、y轴的刻度
hold on
TranMatrix1=[cos(FinalTheta/180*pi) -sin(FinalTheta/180*pi);...
            sin(FinalTheta/180*pi) cos(FinalTheta/180*pi)];
for k=1:length(P(1,:))%坐标转换
    P(:,k)=TranMatrix1*P(:,k);
end
for k=1:length(P(1,:))
    P(1,k)=P(1,k)+FinalVector(1);
    P(2,k)=P(2,k)+FinalVector(2);
end
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4);  
end
for i=1:size(X,2)
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off

%*********************************************************************
%Step4.6：计算MSE,任意一点的真实点为其欧氏距离最近点。

%Step4.6.1：从X点集中寻找Pi点集中每个点的最近点,并得到对应点对
PP=P;
XX=X;
setIn=[PP;XX];                      %两个2×n矩阵(Pi和X)构造一个4×n矩阵矩阵(setIn)，来满足ClosetPointMatch输入参数要求
[setOut1]=ClosetPointMatch1(setIn);  %获取点集Pi的最近点集
Xi=[setOut1(3,:);setOut1(4,:)];
%Step4.6.2：计算XX与PP的误差 
coordinateOffste=Xi-PP;             %计算Xi与Pi的横纵坐标差值
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%计算每个P点和X点之间的欧式距离
end%for j=1:length(coordinateOffste(1,:)) 

MSE=sum(distanceIndiv)/length(distanceIndiv);%P集和X集的平均误差
    

