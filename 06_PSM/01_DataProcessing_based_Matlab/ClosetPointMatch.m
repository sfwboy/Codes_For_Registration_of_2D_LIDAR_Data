

function [set]=ClosetPointMatch(XP)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ�Ѱ�Ҵ�ƥ��㼯ÿ���������㣨�㵽�߶ξ�����С�� 
%���룺1��4��1����XP
%�����1��4��1����set
%���У�XP = |xP|         set=|xP      |
%           |yP|             |yP      |
%           |xX|             |xPCloset|
%           |yX|             |yPCloset|
%���У�(xP,yP)Ϊ��ƥ��㼯P�����ݵ㣬(xX,yX)ģ��㼯X�����ݵ�
%      (xPCloset,yPCloset)λ��ģ��㼯X��Ϊ(xP,yP)һһ��Ӧ�������
%      ���е�����ڵѿ�������ϵ      
%���ߣ�Shaofeng Wu 
%ʱ�䣺2017.11.26
%*********************************************************************
%*********************************************************************
%*********************************************************************


%*********************************************************************
%*********************************************************************
%Step1����ȡ��ʼ����
xP=XP(1,:);         %��ȡP�㼯���е�ĺ�����                     
yP=XP(2,:);         %��ȡP�㼯���е��������   
xX=XP(3,:);         %��ȡX�㼯���е�ĺ�����   
yX=XP(4,:);         %��ȡX�㼯���е��������   
set=[xP;yP;zeros(size(xP));zeros(size(xP))]; %
valueFlag=0;
%*********************************************************************
%*********************************************************************
%Step2����X�㼯��Ѱ�ҵ㼯P�������
for i=1:length(xP)
    rangeMin =90000;
    for j=1:length(xX)-1    %��X�㼯��Ѱ�ҵ㼯P�������
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
        
        if distance<rangeMin    %�����ǰ�����ƥ���ľ���distance��rangeMin ��
            rangeMin=distance;  %�����ΪrangeMin ֵ
            if (valueFlag==1)
                distance1=((xP(i)-xX(j))^2+(yP(i)-yX(j))^2);
                distance2=((xP(i)-xX(j+1))^2+(yP(i)-yX(j+1))^2);
                    if distance1>distance2
                        set(3,i)=xX(j+1);           %���������
                        set(4,i)=yX(j+1);
                    else
                        set(3,i)=xX(j);             %���������
                        set(4,i)=yX(j);
                    end
            elseif (valueFlag==3)
                set(3,i)=xX(j);                     %���������
                set(4,i)=yX(j);
            else
                set(3,i)=xX(j+1);                   %���������
                set(4,i)=yX(j+1);
            end
        end
    end%for j=1:length(xX)   
end%for i=1:length(xP)










