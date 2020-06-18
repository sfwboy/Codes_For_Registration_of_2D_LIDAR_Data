%**************************************************************************
%**************************************************************************
%功能：使用激光雷达产生仿真数据，并且存为.scan文件
%PSM源码中的激光测距仪的数据存储在.scan文件中，数据格式为N×2,如下：
%data（i,1）为激光测距仪的采样角，单位为弧度制
%data（i,2）为激光测距仪在采样角data（i,1）时障碍物的距离，单位为cm
%data中，N=200,0<=data（i,1）<=2*pi
%作者：Shaofeng Wu 
%时间：2019.12.07
%邮箱：shaofeng693@126.com
%**************************************************************************
%**************************************************************************
% function main
clear
clc
%*********************************************************************
%*********************************************************************
%Step1：数据获取和数据预处理
%*********************************************************************
%Step1.1：设置激光扫描仪的前后两个位置，和采样点数
xPre=190;       %激光扫描仪前一时刻位置横坐标
yPre=210;       %激光扫描仪前一时刻位置纵坐标
directionPre=0; %激光扫描仪前一时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
xLate=290;      %激光扫描仪当前时刻位置横坐标
yLate=280;      %激光扫描仪当前时刻位置纵坐标
directionLate=45;%激光扫描仪当前时刻朝向与基准方向的夹角（单位：度），基准方向设为水平向右
numSample=681;  %激光扫描仪360度细分成numSample份
MapName='Map13.bmp';%选择实验地图（场景）
% MapName='Map14.bmp';%选择实验地图（场景）
%*********************************************************************
%Step1.2：获取激光扫描仪前后两个位置的坐标对应的X、P点集
[setOutX]=GetLaserDataPointSet(xPre,yPre,numSample,directionPre,11,MapName);     %初始位置设为（xPre,yPre），返回的数据点集记为模板点集X
[setOutP]=GetLaserDataPointSet(xLate,yLate,numSample,directionLate,22,MapName);%移动位置设为（xLate,yLate），返回的数据点集记为待匹配点集P
%Step1.3：坐标转换，将数据由极坐标系变换为笛卡尔坐标系
X=CoordinateTran2(setOutX);%坐标变换：极坐标变换为笛卡尔坐标
P=CoordinateTran2(setOutP);%坐标变换：极坐标变换为笛卡尔坐标
 %绘制数据点集平面图，可视化激光测距仪采集到的数据
figure(101);
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'rx','MarkerSize',3); 
% plot([P(2,i) P(2,i)],[P(1,i) P(1,i)],'rx','MarkerSize',3);
end
for i=1:size(X,2)
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'gx','MarkerSize',3);  
% plot([X(2,i) X(2,i)],[X(1,i) X(1,i)],'gx','MarkerSize',3); 
end
plot(0,0,'b.','MarkerSize',30 );
hold off


data1=setOutX';
data2=setOutP';

for i=1:size(data1,1)
    dataOut1(i,1)=data1(i,2);
    dataOut1(i,2)=data1(i,1);
    
    dataOut2(i,1)=data2(i,2);
    dataOut2(i,2)=data2(i,1);
end
% 
dlmwrite('Mydata1.scan',dataOut1,'delimiter',' ','precision','%.3f');
dlmwrite('Mydata2.scan',dataOut2,'delimiter',' ','precision','%.3f');

