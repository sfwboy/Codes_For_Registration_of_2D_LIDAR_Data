function [dx,dy,dTheta,sumError]=Registration_of_lines(X,P)
%作者：Shaofeng Wu 
%时间：2019.03.09
%邮箱：shaofeng693@126.com

%*********************************************************************
%*********************************************************************
%Step1：激光测距仪数据的直线段提取
%*********************************************************************

%*********************************************************************
%Step1.1：数据集自适应分割
[SetOut11]=SplitData(X,8,3);%数据点分割
[SetOut22]=SplitData(P,8,3);%数据点分割
%*********************************************************************
%Step1.2：对分割后的数据子集，判断其是否为直线段特征，如果不是则再次分割，然
%         后使用直线特征拟合。
% [SetOut1]=SplitData1(SetOut11,15,3);%对非直线特征数据集只分割一次
% [SetOut2]=SplitData1(SetOut22,15,3);%对非直线特征数据集只分割一次
% [SetOut1]=SplitData2(SetOut11,15,5);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
% [SetOut2]=SplitData2(SetOut22,15,4);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
% [SetOut1]=SplitData2(SetOut11,15,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
% [SetOut2]=SplitData2(SetOut22,15,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
[SetOut1]=SplitData2(SetOut11,10,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集
[SetOut2]=SplitData2(SetOut22,10,3);%对非直线特征数据集分割n次，直到所有子集都为直线特征数据集

% 显示数据集分割和处理后的情况（主要是调试时用）
% figure(102);
% hold on
% for j=1:size(SetOut1,3)
%     for i=1:size(SetOut1,2)
%         if SetOut1(1,i,j)~=0 || SetOut1(2,i,j)~=0
%             plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'gx','MarkerSize',3);    
%         end
%     end
% end
% for j=1:size(SetOut2,3)
%     for i=1:size(SetOut2,2)
%         if SetOut2(1,i,j)~=0 || SetOut2(2,i,j)~=0
%             plot([SetOut2(1,i,j) SetOut2(1,i,j)],[SetOut2(2,i,j) SetOut2(2,i,j)],'rx','MarkerSize',3);    
%         end
%     end
% end
% hold off




%*********************************************************************
%Step1.3：直线段拟合
% [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1);%利用最小二乘法拟合直线，效果很好
% [Result_Ls2,NumFlag2]=LeastSquareLine1(SetOut2);%利用最小二乘法拟合直线，效果很好
% [Result_Ls1,NumFlag1]=LineFittingRANSAC(SetOut1);%利用RANSAC法拟合直线，效果一般
% [Result_Ls2,NumFlag2]=LineFittingRANSAC(SetOut2);%利用RANSAC法拟合直线，效果一般
[Result_Ls1,NumFlag1]=LineFitting(SetOut1);%RANSAC法改进为遍历数据集，拟合直线，效果很好
[Result_Ls2,NumFlag2]=LineFitting(SetOut2);%RANSAC法改进为遍历数据集，拟合直线，效果很好

%显示最终的拟合情况
% LinePoint1=zeros(2,2,size(SetOut1,3));%收集直线段的端点
% LinePoint2=zeros(2,2,size(SetOut2,3));%收集直线段的端点
% figure(103);
% hold on
% for j=1:size(SetOut1,3)
%     for i=1:size(SetOut1,2)
%         if SetOut1(1,i,j)~=0 || SetOut1(2,i,j)~=0
%             plot([SetOut1(1,i,j) SetOut1(1,i,j)],[SetOut1(2,i,j) SetOut1(2,i,j)],'gx','MarkerSize',3);    
%         end
%     end
% end
% for j=1:size(SetOut2,3)
%     for i=1:size(SetOut2,2)
%         if SetOut2(1,i,j)~=0 || SetOut2(2,i,j)~=0
%             plot([SetOut2(1,i,j) SetOut2(1,i,j)],[SetOut2(2,i,j) SetOut2(2,i,j)],'rx','MarkerSize',3);    
%         end
%     end
% end

for j=1:size(SetOut1,3)
    if abs(SetOut1(1,1,j)-SetOut1(1,NumFlag1(j),j))>5%线段首尾点横坐标距离作为阈值，保证获取端点贴合数据
%         plot([SetOut1(1,1,j) SetOut1(1,NumFlag1(j),j)],...
%         [Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j) ...
%         Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j)],'b-x','LineWidth',1,'MarkerSize',3);
        LinePoint1(1,1,j)=SetOut1(1,1,j);%起始点横坐标
        LinePoint1(1,2,j)=SetOut1(1,NumFlag1(j),j);%结束点横坐标
        LinePoint1(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%起始点纵坐标
        LinePoint1(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%结束点纵坐标
    else
%         plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
%         (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
%         [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-x','LineWidth',1,'MarkerSize',3);
        LinePoint1(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%起始点横坐标
        LinePoint1(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%结束点横坐标
        LinePoint1(2,1,j)=SetOut1(2,1,j);%起始点纵坐标
        LinePoint1(2,2,j)=SetOut1(2,NumFlag1(j),j);%结束点纵坐标
    end 
end
for j=1:size(SetOut2,3)
    if abs(SetOut2(1,1,j)-SetOut2(1,NumFlag2(j),j))>5
%         plot([SetOut2(1,1,j) SetOut2(1,NumFlag2(j),j)],...l
%         [Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j) ...
%         Result_Ls2(1,j)*SetOut2(1,NumFlag2(j),j)+Result_Ls2(2,j)],'k-x','LineWidth',1,'MarkerSize',3);
        LinePoint2(1,1,j)=SetOut2(1,1,j);%起始点横坐标
        LinePoint2(1,2,j)=SetOut2(1,NumFlag2(j),j);%结束点横坐标
        LinePoint2(2,1,j)=Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j);%起始点纵坐标
        LinePoint2(2,2,j)=Result_Ls2(1,j)*SetOut2(1,NumFlag2(j),j)+Result_Ls2(2,j);%结束点纵坐标
    else
%         plot([(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j) ...
%         (SetOut2(2,NumFlag2(j),j)-Result_Ls2(2,j))/Result_Ls2(1,j)],...
%         [SetOut2(2,1,j) SetOut2(2,NumFlag2(j),j)],'k-x','LineWidth',1,'MarkerSize',3);
        LinePoint2(1,1,j)=(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j);%起始点横坐标
        LinePoint2(1,2,j)=(SetOut2(2,NumFlag2(j),j)-Result_Ls2(2,j))/Result_Ls2(1,j);%结束点横坐标
        LinePoint2(2,1,j)=SetOut2(2,1,j);%起始点纵坐标
%         LinePoint2(2,2,j)=SetOut2(2,NumFlag1(j),j);%结束点纵坐标
        LinePoint2(2,2,j)=SetOut2(2,NumFlag2(j),j);%结束点纵坐标
    end 
end
% hold off

%使用直线的端点显示直线
% figure(104);
% hold on
% for i=1:size(LinePoint1,3)
%    plot([LinePoint1(1,1,i) LinePoint1(1,2,i)],[LinePoint1(2,1,i) LinePoint1(2,2,i)],'g-','LineWidth',2);%绘制直线
%    %绘制端点和原点的连线
%    plot([LinePoint1(1,1,i) 0],[LinePoint1(2,1,i) 0],'b--');%绘制直线  
%    plot([0 LinePoint1(1,2,i)],[0 LinePoint1(2,2,i)],'y--');%绘制直线   
% end
% for i=1:size(LinePoint2,3)
%    plot([LinePoint2(1,1,i) LinePoint2(1,2,i)],[LinePoint2(2,1,i) LinePoint2(2,2,i)],'r-','LineWidth',2);%绘制直线
%    %绘制端点和原点的连线
%    plot([LinePoint2(1,1,i) 0],[LinePoint2(2,1,i) 0],'r--');%绘制直线  
%    plot([0 LinePoint2(1,2,i)],[0 LinePoint2(2,2,i)],'m--');%绘制直线
% end
% hold off


[thetaD1]=LineDomain(LinePoint1,Result_Ls1);
%保持同帧数据的不同直线定义域互斥
if size(thetaD1,2)>2
    for i=1:size(thetaD1,2)-1
%         if (thetaD1(1,i)<=thetaD1(1,i+1)<=thetaD1(2,i))||(thetaD1(1,i+1)<=thetaD1(1,i)<=thetaD1(2,i+1))
        if (thetaD1(1,i)<thetaD1(1,i+1))&&(thetaD1(1,i+1)<thetaD1(2,i))&&(thetaD1(2,i)<thetaD1(2,i+1))
            thetaD1(2,i)=thetaD1(1,i+1)-0.001;%0.001/pi*180=0.057度
        end
    end
end
[thetaD2]=LineDomain(LinePoint2,Result_Ls2);
%保持同帧数据的不同直线定义域互斥
if size(thetaD2,2)>2
    for i=1:size(thetaD2,2)-1
%         if (thetaD1(1,i)<=thetaD1(1,i+1)<=thetaD1(2,i))||(thetaD1(1,i+1)<=thetaD1(1,i)<=thetaD1(2,i+1))
        if (thetaD2(1,i)<thetaD2(1,i+1))&&(thetaD2(1,i+1)<thetaD2(2,i))&&(thetaD2(2,i)<thetaD2(2,i+1))
            thetaD2(2,i)=thetaD2(1,i+1)-0.001;%0.001/pi*180=0.057度
        end
    end
end


num=0;
for i=1:size(thetaD1,2)
    thetaA1=thetaD1(1,i);
    thetaB1=thetaD1(2,i);
    for j=1:size(thetaD2,2)
        thetaA2=thetaD2(1,j);
        thetaB2=thetaD2(2,j);
        if((thetaA2<=thetaA1)&&(thetaA1<=thetaB2 && thetaB2<=thetaB1))||...
                ((thetaA1<=thetaA2)&&(thetaA2<=thetaB1 && thetaB1<=thetaB2))||...
                (thetaA2<=thetaA1 && thetaB1<=thetaB2)||...
                (thetaA1<=thetaA2 && thetaB2<=thetaB1)
            num=num+1;
            reDomain(1,num)=max(thetaA1,thetaA2);%重合直线的起始角度φai
            reDomain(2,num)=min(thetaB1,thetaB2);%重合直线的结束角度φbi
            reDomain(3,num)=thetaD1(3,i);        %直线1法向量与x轴正方向的夹角αi
            reDomain(4,num)=thetaD1(4,i);        %原点到直线1的距离li
            reDomain(5,num)=thetaD2(3,j);        %直线2法向量与x轴正方向的夹角αi'
            reDomain(6,num)=thetaD2(4,j);        %原点到直线2的距离li'
            reDomain(7,num)=thetaD1(5,i);        %直线1的参数a,(y=ax+b)
            reDomain(8,num)=thetaD1(6,i);        %直线1的参数b,(y=ax+b)
            reDomain(9,num)=thetaD2(5,j);        %直线2的参数a,(y=ax+b)
            reDomain(10,num)=thetaD2(6,j);       %直线2,的参数b,(y=ax+b)
        end
    end
    
end
% figure(105)
% hold on
% for i=1:size(LinePoint1,3)
%    plot([LinePoint1(1,1,i) LinePoint1(1,2,i)],[LinePoint1(2,1,i) LinePoint1(2,2,i)],'g-','LineWidth',2);%绘制直线    
% end
% for i=1:size(LinePoint2,3)
%    plot([LinePoint2(1,1,i) LinePoint2(1,2,i)],[LinePoint2(2,1,i) LinePoint2(2,2,i)],'r-','LineWidth',2);%绘制直线
% end
%     for i=1:size(reDomain,2)
%         k1=tan(reDomain(1,i));%经过原点，经过直线定义域起点的直线斜率
%         k2=tan(reDomain(2,i));%经过原点，经过直线定义域结束点的直线斜率
%         %求解直线1与上述两条直线的交点
%         %l1:y=ax+b,k1:y=k1x,k2:y=k2x
%         x11=reDomain(8,i)/(k1-reDomain(7,i));
%         y11=x11*k1;
%         x12=reDomain(8,i)/(k2-reDomain(7,i));
%         y12=x12*k2;
%         
%         %求解直线2与上述两条直线的交点
%         %l1:y=ax+b,k1:y=k1x,k2:y=k2x
%         x21=reDomain(10,i)/(k1-reDomain(9,i));
%         y21=x21*k1;
%         x22=reDomain(10,i)/(k2-reDomain(9,i));
%         y22=x22*k2;
%         
%         plot([0 x11],[0 y11],'m--');%绘制直线
%         plot([0 x12],[0 y12],'k--');%绘制直线
%         plot([0 x21],[0 y21],'m--');%绘制直线
%         plot([0 x22],[0 y22],'k--');%绘制直线
%         
%         plot([x11 x21],[y11 y21],'m-','LineWidth',2);%绘制直线
%         plot([x12 x22],[y12 y22],'k-','LineWidth',2);%绘制直线
%     end
% hold off


sumError=0;

%总的面积误差
factor=0;
for i=1:size(reDomain,2)
    factor=factor+reDomain(2,i)-reDomain(1,i);
    if reDomain(3,i)==reDomain(5,i)
        tempErrorAi=reDomain(4,i)^2*tan(reDomain(1,i)-reDomain(3,i))+...
                    reDomain(6,i)^2*tan(reDomain(1,i)-reDomain(5,i))+...
                    2*reDomain(4,i)*reDomain(6,i)*tan(reDomain(1,i)-reDomain(3,i));
        tempErrorBi=reDomain(4,i)^2*tan(reDomain(2,i)-reDomain(3,i))+...
                    reDomain(6,i)^2*tan(reDomain(2,i)-reDomain(5,i))+...
                    2*reDomain(4,i)*reDomain(6,i)*tan(reDomain(2,i)-reDomain(3,i));
        tempError=tempErrorBi-tempErrorAi;
        sumError=sumError+tempError;
    else
        tempErrorAi=reDomain(4,i)^2*tan(reDomain(1,i)-reDomain(3,i))+...
                    reDomain(6,i)^2*tan(reDomain(1,i)-reDomain(5,i))+...
                    2*reDomain(4,i)*reDomain(6,i)*...
                    log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
                    sin(reDomain(3,i)-reDomain(5,i));
        tempErrorBi=reDomain(4,i)^2*tan(reDomain(2,i)-reDomain(3,i))+...
                    reDomain(6,i)^2*tan(reDomain(2,i)-reDomain(5,i))+...
                    2*reDomain(4,i)*reDomain(6,i)*...
                    log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
                    sin(reDomain(3,i)-reDomain(5,i));
        tempError=tempErrorBi-tempErrorAi;
        sumError=sumError+tempError; 
%         tempError=abs(tempErrorBi)-abs(tempErrorAi);
%         sumError=sumError+abs(tempError); 
    end
end
%添加系数
factor=2*pi/(factor^2);
sumError=sumError*factor;

%求导数
dx=0;
dy=0;
dTheta=0;
for i=1:size(reDomain,2)
    if reDomain(3,i)==reDomain(5,i)
        %在求解dE/dθ时，reDomain(3,i)==reDomain(5,i)会导致dE/dθ区域无穷大，导致无法往下迭代，所以
        %人为加一个偏差，使得两者不在相等
        reDomain(5,i)=reDomain(5,i)+1/180*pi;
        %求关于x（横坐标）的导数
        dxAi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(1,i)-reDomain(3,i)));
        dxBi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(2,i)-reDomain(3,i)));
        dxTemp=dxBi-dxAi;
        dx=dx+dxTemp;
        %求关于y（纵坐标）的导数
        dyAi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(1,i)-reDomain(3,i)));
        dyBi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(2,i)-reDomain(3,i)));
        dyTemp=dyBi-dyAi;
        dy=dy+dyTemp;
        %求关于theta（旋转角度）的导数
        dThetaAi=-reDomain(6,i)^2/(cos(reDomain(1,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(1,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*tan(reDomain(1,i)-reDomain(3,i)));
        dThetaBi=-reDomain(6,i)^2/(cos(reDomain(2,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(2,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*tan(reDomain(2,i)-reDomain(3,i)));
        dThetaTemp=dThetaBi-dThetaAi;
        dTheta=dTheta+dThetaTemp;
    else
        %求关于x（横坐标）的导数
        dxAi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dxBi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dxTemp=dxBi-dxAi;
        dx=dx+dxTemp;
        %求关于y（纵坐标）的导数
        dyAi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dyBi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dyTemp=dyBi-dyAi;
        dy=dy+dyTemp;
        %求关于theta（旋转角度）的导数
        dThetaAi=-reDomain(6,i)^2/(cos(reDomain(1,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(1,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*...
                 log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
                 sin(reDomain(3,i)-reDomain(5,i)));
        dThetaBi=-reDomain(6,i)^2/(cos(reDomain(2,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(2,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*...
                 log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
                 sin(reDomain(3,i)-reDomain(5,i)));
        dThetaTemp=dThetaBi-dThetaAi;
        dTheta=dTheta+dThetaTemp;     
    end
end

% dx=dx*LearnningRate;
% dy=dy*LearnningRate;
% dTheta=dTheta*LearnningRate;


end