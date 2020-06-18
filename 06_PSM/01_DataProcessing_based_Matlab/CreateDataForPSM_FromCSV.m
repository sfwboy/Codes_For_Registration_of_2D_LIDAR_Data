clc
clear
%**************************************************************************
%**************************************************************************
%功能：将csv文件数据转化为.scan文件数据。
%PSM源码中的激光测距仪的数据存储在.scan文件中，数据格式为N×2,如下：
%data（i,1）为激光测距仪的采样角，单位为弧度制
%data（i,2）为激光测距仪在采样角data（i,1）时障碍物的距离，单位为cm
%data中，N=200,0<=data（i,1）<=2*pi
%作者：Shaofeng Wu 
%时间：2019.12.07
%邮箱：shaofeng693@126.com
%**************************************************************************
%**************************************************************************
%*********************************************************************
%Step1.1：从存储的csv文件中读取数据
setOutX = csvread('Exp1_X_NoiseAdd_5dB.csv');
setOutP = csvread('Exp1_P_NoiseAdd_5dB.csv');
%*********************************************************************
%Step1.2：坐标转换，将数据由极坐标系变换为笛卡尔坐标系
X=CoordinateTran2(setOutX);%坐标变换：极坐标变换为笛卡尔坐标
P=CoordinateTran2(setOutP);%坐标变换：极坐标变换为笛卡尔坐标
%*********************************************************************
%Step1.3：坐标转换，将数据由极坐标系变换为笛卡尔坐标系
figure(101);
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
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
dlmwrite('Exp1_X_NoiseAdd_5dB.scan',dataOut1,'delimiter',' ','precision','%.3f');
dlmwrite('Exp1_P_NoiseAdd_5dB.scan',dataOut2,'delimiter',' ','precision','%.3f');