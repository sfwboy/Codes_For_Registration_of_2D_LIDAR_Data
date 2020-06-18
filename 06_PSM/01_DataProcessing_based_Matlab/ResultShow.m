clc
clear
%*********************************************************************
%Step1.1：从存储的csv文件中读取数据
setOutX = csvread('Exp4_X_NoiseAdd_5dB.csv');
setOutP = csvread('Exp4_P_NoiseAdd_5dB.csv');
FinalTheta=-26.6738/180*pi;               %根据PSM源码得到的旋转角度
FinalVector=[24.0180,-35.7921]; %根据PSM源码得到的平移矩阵
%*********************************************************************
%Step1.2：坐标转换，将数据由极坐标系变换为笛卡尔坐标系
X=CoordinateTran(setOutX);%坐标变换：极坐标变换为笛卡尔坐标
P=CoordinateTran(setOutP);%坐标变换：极坐标变换为笛卡尔坐标 
%Step1.3：显示原始数据
figure(1111111);
% title('原始激光数据点');    %标题
% xlabel('水平距离(cm)');     %x轴
% ylabel('竖直距离(cm)');     %y轴
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4); 
     plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off

% FinalVector=[20,-20];
TranMatrix1=[cos(FinalTheta) -sin(FinalTheta);...
            sin(FinalTheta) cos(FinalTheta)];
for k=1:length(P(1,:))%坐标转换
    P1(:,k)=TranMatrix1*P(:,k);
end
for k=1:length(P1(1,:))
    P1(1,k)=P1(1,k)+FinalVector(1);
    P1(2,k)=P1(2,k)+FinalVector(2);
end

figure(1111112);
% title('原始激光数据点');    %标题
% xlabel('水平距离(cm)');     %x轴
% ylabel('竖直距离(cm)');     %y轴
hold on
for i=1:size(P1,2)
    plot([P1(1,i) P1(1,i)],[P1(2,i) P1(2,i)],'r.','MarkerSize',4); 
     plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off








% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%求解MSE,任意一点的真实点为其欧氏距离最近点。
% %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.2.1：从X点集中寻找Pi点集中每个点的最近点,并得到对应点对
PP=P1;
XX=X;
setIn=[PP;XX];                       %两个2×n矩阵(Pi和X)构造一个4×n矩阵矩阵(setIn)，来满足ClosetPointMatch输入参数要求
[setOut1]=ClosetPointMatch1(setIn);  %获取点集Pi的最近点集
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.2.2：计算X与Pi的误差 
coordinateOffste=Xi-PP;             %计算Xi与Pi的横纵坐标差值
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%计算每个P点和X点之间的欧式距离
end%for j=1:length(coordinateOffste(1,:)) 

MSE=sum(distanceIndiv)/length(distanceIndiv);%P集和X集的平均误差