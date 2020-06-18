function [thetaD1]=LineDomain(LinePoint1,Result_Ls1)
%作者：Shaofeng Wu 
%时间：2018.10.19
%邮箱：shaofeng693@126.com
    for i=1:size(LinePoint1,3)
    %输出为thetaD1，
    %thetaD1(1,i)=直线定义域起点
    %thetaD1(2,i)=直线定义域结束点
    %thetaD1(3,i)=直线法向量与x轴正方向的夹角
    %thetaD1(4,i)=原点到直线的距离
    %thetaD1(5,i)=直线的参数a（y=ax+b）
    %thetaD1(6,i)=直线的参数a（y=ax+b）
    %***********************************************************************
    %***********************************************************************
    %Step1:计算直线段组端点和x轴正方向的夹角（即直线段夹角的范围）
    %用向量（1,0）代表x轴正方向,使用余弦定理求解两个向量的夹角
    %起始点
    if (LinePoint1(1,1,i)<0 && (LinePoint1(2,1,i)<0))  %端点位于第三象限
        thetaD1(1,i)=2*pi-acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    elseif(LinePoint1(1,1,i)>0 && (LinePoint1(2,1,i)<0))  %端点位于第四象限
        thetaD1(1,i)=2*pi-acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    else   %端点位于第一、二象限
        thetaD1(1,i)=acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    end
    %结束点
    if (LinePoint1(1,2,i)<0 && (LinePoint1(2,2,i)<0))  %端点位于第三象限
        thetaD1(2,i)=2*pi-acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    elseif(LinePoint1(1,2,i)>0 && (LinePoint1(2,2,i)<0))  %端点位于第四象限
        thetaD1(2,i)=2*pi-acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    else   %端点位于第一、二象限
        thetaD1(2,i)=acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    end
   
    %***********************************************************************
    %***********************************************************************
    %Step2:计算直线法向量与x轴正方向的夹角，以及原点到直线的距离

    %%%计算直线法向量与x轴正方向的夹角
    %用向量（1,0）代表x轴正方向,使用余弦定理求解两个向量的夹角
    dirVecL1(:,i)=LinePoint1(:,2,i)-LinePoint1(:,1,i);%线段起始点相减即为线段的方向向量
    norVecL1(1,i)=1;
    norVecL1(2,i)=-dirVecL1(1,i)/dirVecL1(2,i);%直线段法向量
    %plot([0 norVecL1(1,i)],[0 norVecL1(2,i)],'y-','LineWidth',2);
    %thetaD1(3,i)=acos(norVecL1(1,i)/(sqrt(1+norVecL1(2,i)^2)));%直线法向量与x轴正方向的夹角
    X_xy=[1;0];
    tempTheta=acos(dot(X_xy,norVecL1(:,i))/(norm(X_xy)*norm(norVecL1(:,i))));%直线法向量与x轴正方向的夹角
    if(norVecL1(2,i)<0)
        thetaD1(3,i)=2*pi-tempTheta;
    else
        thetaD1(3,i)=tempTheta;
    end
   
    if ~(((thetaD1(3,i)-pi/2)<=thetaD1(1,i))&&((thetaD1(3,i)+pi/2)>=thetaD1(2,i)))%如果不满足条件，则旋转180度
        if ((thetaD1(3,i)-2*pi-pi/2)<=thetaD1(1,i))&&((thetaD1(3,i)-2*pi+pi/2)>=thetaD1(2,i))%旋转前判断是否该角度为x轴附近
            thetaD1(3,i)=thetaD1(3,i);%保持不变
        else
            thetaD1(3,i)=thetaD1(3,i)+pi;%旋转180度
        end
    end
    if(thetaD1(3,i)>2*pi)
        thetaD1(3,i)=thetaD1(3,i)-2*pi;%减掉360，只是为了方便测试，理论上减不减都一样,因为sin（phi+-2kpi）=sin(phi),cos、tan、cot等同。
    end
    %***********************************************************************
    %***********************************************************************
    %Step3:计算原点到直线的距离，使用点到直线距离公式
    %直线公式为：Ax+By+C=0,点坐标为(x0,y0),则
    %d=|(A*x0+B*y0+C)/(sqrt(A^2+B^2)|
    %该程序中，直线公式为y=ax+b -->ax-y+b=0,
    %该程序中，点为原点，坐标为(0,0)
    %所以A=a,B=-1,C=b;d=|b/sqrt(a^2+1)|
    thetaD1(4,i)=abs(Result_Ls1(2,i)/sqrt(Result_Ls1(1,i)^2+1));%原点到直线的距离
    thetaD1(5,i)=Result_Ls1(1,i);%直线的参数a
    thetaD1(6,i)=Result_Ls1(2,i);%直线参数b
    end
    for i=1:size(thetaD1,2)%拆分起始点角度大于结束点角度的直线
       if(thetaD1(1,i)>thetaD1(2,i)) 
          thetaD1(2,i)=2*pi+thetaD1(2,i);

       end
    end
end
