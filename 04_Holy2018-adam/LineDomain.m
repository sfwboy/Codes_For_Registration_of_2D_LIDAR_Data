function [thetaD1]=LineDomain(LinePoint1,Result_Ls1)
%���ߣ�Shaofeng Wu 
%ʱ�䣺2018.10.19
%���䣺shaofeng693@126.com
    for i=1:size(LinePoint1,3)
    %���ΪthetaD1��
    %thetaD1(1,i)=ֱ�߶��������
    %thetaD1(2,i)=ֱ�߶����������
    %thetaD1(3,i)=ֱ�߷�������x��������ļн�
    %thetaD1(4,i)=ԭ�㵽ֱ�ߵľ���
    %thetaD1(5,i)=ֱ�ߵĲ���a��y=ax+b��
    %thetaD1(6,i)=ֱ�ߵĲ���a��y=ax+b��
    %***********************************************************************
    %***********************************************************************
    %Step1:����ֱ�߶���˵��x��������ļнǣ���ֱ�߶μнǵķ�Χ��
    %��������1,0������x��������,ʹ�����Ҷ���������������ļн�
    %��ʼ��
    if (LinePoint1(1,1,i)<0 && (LinePoint1(2,1,i)<0))  %�˵�λ�ڵ�������
        thetaD1(1,i)=2*pi-acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    elseif(LinePoint1(1,1,i)>0 && (LinePoint1(2,1,i)<0))  %�˵�λ�ڵ�������
        thetaD1(1,i)=2*pi-acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    else   %�˵�λ�ڵ�һ��������
        thetaD1(1,i)=acos((LinePoint1(1,1,i)*1+LinePoint1(2,1,i)*0)/...
                        sqrt(LinePoint1(1,1,i)^2+LinePoint1(2,1,i)^2)*sqrt(1^2+0^2));
    end
    %������
    if (LinePoint1(1,2,i)<0 && (LinePoint1(2,2,i)<0))  %�˵�λ�ڵ�������
        thetaD1(2,i)=2*pi-acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    elseif(LinePoint1(1,2,i)>0 && (LinePoint1(2,2,i)<0))  %�˵�λ�ڵ�������
        thetaD1(2,i)=2*pi-acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    else   %�˵�λ�ڵ�һ��������
        thetaD1(2,i)=acos((LinePoint1(1,2,i)*1+LinePoint1(2,2,i)*0)/...
                        sqrt(LinePoint1(1,2,i)^2+LinePoint1(2,2,i)^2)*sqrt(1^2+0^2));
    end
   
    %***********************************************************************
    %***********************************************************************
    %Step2:����ֱ�߷�������x��������ļнǣ��Լ�ԭ�㵽ֱ�ߵľ���

    %%%����ֱ�߷�������x��������ļн�
    %��������1,0������x��������,ʹ�����Ҷ���������������ļн�
    dirVecL1(:,i)=LinePoint1(:,2,i)-LinePoint1(:,1,i);%�߶���ʼ�������Ϊ�߶εķ�������
    norVecL1(1,i)=1;
    norVecL1(2,i)=-dirVecL1(1,i)/dirVecL1(2,i);%ֱ�߶η�����
    %plot([0 norVecL1(1,i)],[0 norVecL1(2,i)],'y-','LineWidth',2);
    %thetaD1(3,i)=acos(norVecL1(1,i)/(sqrt(1+norVecL1(2,i)^2)));%ֱ�߷�������x��������ļн�
    X_xy=[1;0];
    tempTheta=acos(dot(X_xy,norVecL1(:,i))/(norm(X_xy)*norm(norVecL1(:,i))));%ֱ�߷�������x��������ļн�
    if(norVecL1(2,i)<0)
        thetaD1(3,i)=2*pi-tempTheta;
    else
        thetaD1(3,i)=tempTheta;
    end
   
    if ~(((thetaD1(3,i)-pi/2)<=thetaD1(1,i))&&((thetaD1(3,i)+pi/2)>=thetaD1(2,i)))%�������������������ת180��
        if ((thetaD1(3,i)-2*pi-pi/2)<=thetaD1(1,i))&&((thetaD1(3,i)-2*pi+pi/2)>=thetaD1(2,i))%��תǰ�ж��Ƿ�ýǶ�Ϊx�ḽ��
            thetaD1(3,i)=thetaD1(3,i);%���ֲ���
        else
            thetaD1(3,i)=thetaD1(3,i)+pi;%��ת180��
        end
    end
    if(thetaD1(3,i)>2*pi)
        thetaD1(3,i)=thetaD1(3,i)-2*pi;%����360��ֻ��Ϊ�˷�����ԣ������ϼ�������һ��,��Ϊsin��phi+-2kpi��=sin(phi),cos��tan��cot��ͬ��
    end
    %***********************************************************************
    %***********************************************************************
    %Step3:����ԭ�㵽ֱ�ߵľ��룬ʹ�õ㵽ֱ�߾��빫ʽ
    %ֱ�߹�ʽΪ��Ax+By+C=0,������Ϊ(x0,y0),��
    %d=|(A*x0+B*y0+C)/(sqrt(A^2+B^2)|
    %�ó����У�ֱ�߹�ʽΪy=ax+b -->ax-y+b=0,
    %�ó����У���Ϊԭ�㣬����Ϊ(0,0)
    %����A=a,B=-1,C=b;d=|b/sqrt(a^2+1)|
    thetaD1(4,i)=abs(Result_Ls1(2,i)/sqrt(Result_Ls1(1,i)^2+1));%ԭ�㵽ֱ�ߵľ���
    thetaD1(5,i)=Result_Ls1(1,i);%ֱ�ߵĲ���a
    thetaD1(6,i)=Result_Ls1(2,i);%ֱ�߲���b
    end
    for i=1:size(thetaD1,2)%�����ʼ��Ƕȴ��ڽ�����Ƕȵ�ֱ��
       if(thetaD1(1,i)>thetaD1(2,i)) 
          thetaD1(2,i)=2*pi+thetaD1(2,i);

       end
    end
end
