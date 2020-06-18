

function [set]=ClosetPointMatch(XP)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：寻找待匹配点集每个点的最近点（点到线段距离最小） 
%输入：1个4×1矩阵XP
%输出：1个4×1矩阵set
%其中，XP = |xP|         set=|xP      |
%           |yP|             |yP      |
%           |xX|             |xPCloset|
%           |yX|             |yPCloset|
%其中，(xP,yP)为待匹配点集P的数据点，(xX,yX)模板点集X的数据点
%      (xPCloset,yPCloset)位于模板点集X，为(xP,yP)一一对应的最近点
%      所有点均处于笛卡尔坐标系      
%作者：Shaofeng Wu 
%时间：2017.11.26
%*********************************************************************
%*********************************************************************
%*********************************************************************


%*********************************************************************
%*********************************************************************
%Step1：获取初始数据
xP=XP(1,:);         %获取P点集所有点的横坐标                     
yP=XP(2,:);         %获取P点集所有点的纵坐标   
xX=XP(3,:);         %获取X点集所有点的横坐标   
yX=XP(4,:);         %获取X点集所有点的纵坐标   
set=[xP;yP;zeros(size(xP));zeros(size(xP))]; %
valueFlag=0;
%*********************************************************************
%*********************************************************************
%Step2：在X点集中寻找点集P的最近点
for i=1:length(xP)
    rangeMin =90000;
    for j=1:length(xX)-1    %在X点集中寻找点集P的最近点
        xAP=xX(j)-xP(i);
        yAP=yX(j)-yP(i);
        xAB=xX(j)-xX(j+1);
        yAB=yX(j)-yX(j+1);
        r_Flag=(xAB*xAP+yAB*yAP)/(xAB^2+yAB^2);
        if (r_Flag>0)&&(r_Flag<1)
            A=[(yX(j)-yX(j+1))/(xX(j)-xX(j+1))];
            B=-1;
            C=yX(j+1)-xX(j+1)*A;
            distance=abs(A*xP(i)+B*yP(i)+C)/sqrt(A^2+B^2);
            valueFlag=1;
        elseif (r_Flag>1)||(r_Flag==1)
           distance=sqrt((xP(i)-xX(j+1))^2+(yP(i)-yX(j+1))^2);
           valueFlag=2;
        else
            distance=sqrt((xP(i)-xX(j))^2+(yP(i)-yX(j))^2);
            valueFlag=3;
        end
        
        if distance<rangeMin    %如果当前点的与匹配点的距离distance比rangeMin ，
            rangeMin=distance;  %则更新为rangeMin 值
            if (valueFlag==1)
                distance1=((xP(i)-xX(j))^2+(yP(i)-yX(j))^2);
                distance2=((xP(i)-xX(j+1))^2+(yP(i)-yX(j+1))^2);
                    if distance1>distance2
                        set(3,i)=xX(j+1);           %更新最近点
                        set(4,i)=yX(j+1);
                    else
                        set(3,i)=xX(j);             %更新最近点
                        set(4,i)=yX(j);
                    end
            elseif (valueFlag==3)
                set(3,i)=xX(j);                     %更新最近点
                set(4,i)=yX(j);
            else
                set(3,i)=xX(j+1);                   %更新最近点
                set(4,i)=yX(j+1);
            end
        end
    end%for j=1:length(xX)   
end%for i=1:length(xP)










