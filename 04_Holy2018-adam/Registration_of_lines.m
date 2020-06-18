function [dx,dy,dTheta,sumError]=Registration_of_lines(X,P)
%���ߣ�Shaofeng Wu 
%ʱ�䣺2019.03.09
%���䣺shaofeng693@126.com

%*********************************************************************
%*********************************************************************
%Step1�������������ݵ�ֱ�߶���ȡ
%*********************************************************************

%*********************************************************************
%Step1.1�����ݼ�����Ӧ�ָ�
[SetOut11]=SplitData(X,8,3);%���ݵ�ָ�
[SetOut22]=SplitData(P,8,3);%���ݵ�ָ�
%*********************************************************************
%Step1.2���Էָ��������Ӽ����ж����Ƿ�Ϊֱ�߶�����������������ٴηָȻ
%         ��ʹ��ֱ��������ϡ�
% [SetOut1]=SplitData1(SetOut11,15,3);%�Է�ֱ���������ݼ�ֻ�ָ�һ��
% [SetOut2]=SplitData1(SetOut22,15,3);%�Է�ֱ���������ݼ�ֻ�ָ�һ��
% [SetOut1]=SplitData2(SetOut11,15,5);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
% [SetOut2]=SplitData2(SetOut22,15,4);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
% [SetOut1]=SplitData2(SetOut11,15,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
% [SetOut2]=SplitData2(SetOut22,15,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
[SetOut1]=SplitData2(SetOut11,10,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
[SetOut2]=SplitData2(SetOut22,10,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�

% ��ʾ���ݼ��ָ�ʹ������������Ҫ�ǵ���ʱ�ã�
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
%Step1.3��ֱ�߶����
% [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1);%������С���˷����ֱ�ߣ�Ч���ܺ�
% [Result_Ls2,NumFlag2]=LeastSquareLine1(SetOut2);%������С���˷����ֱ�ߣ�Ч���ܺ�
% [Result_Ls1,NumFlag1]=LineFittingRANSAC(SetOut1);%����RANSAC�����ֱ�ߣ�Ч��һ��
% [Result_Ls2,NumFlag2]=LineFittingRANSAC(SetOut2);%����RANSAC�����ֱ�ߣ�Ч��һ��
[Result_Ls1,NumFlag1]=LineFitting(SetOut1);%RANSAC���Ľ�Ϊ�������ݼ������ֱ�ߣ�Ч���ܺ�
[Result_Ls2,NumFlag2]=LineFitting(SetOut2);%RANSAC���Ľ�Ϊ�������ݼ������ֱ�ߣ�Ч���ܺ�

%��ʾ���յ�������
% LinePoint1=zeros(2,2,size(SetOut1,3));%�ռ�ֱ�߶εĶ˵�
% LinePoint2=zeros(2,2,size(SetOut2,3));%�ռ�ֱ�߶εĶ˵�
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
    if abs(SetOut1(1,1,j)-SetOut1(1,NumFlag1(j),j))>5%�߶���β������������Ϊ��ֵ����֤��ȡ�˵���������
%         plot([SetOut1(1,1,j) SetOut1(1,NumFlag1(j),j)],...
%         [Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j) ...
%         Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j)],'b-x','LineWidth',1,'MarkerSize',3);
        LinePoint1(1,1,j)=SetOut1(1,1,j);%��ʼ�������
        LinePoint1(1,2,j)=SetOut1(1,NumFlag1(j),j);%�����������
        LinePoint1(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%��ʼ��������
        LinePoint1(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%������������
    else
%         plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
%         (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
%         [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-x','LineWidth',1,'MarkerSize',3);
        LinePoint1(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%��ʼ�������
        LinePoint1(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%�����������
        LinePoint1(2,1,j)=SetOut1(2,1,j);%��ʼ��������
        LinePoint1(2,2,j)=SetOut1(2,NumFlag1(j),j);%������������
    end 
end
for j=1:size(SetOut2,3)
    if abs(SetOut2(1,1,j)-SetOut2(1,NumFlag2(j),j))>5
%         plot([SetOut2(1,1,j) SetOut2(1,NumFlag2(j),j)],...l
%         [Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j) ...
%         Result_Ls2(1,j)*SetOut2(1,NumFlag2(j),j)+Result_Ls2(2,j)],'k-x','LineWidth',1,'MarkerSize',3);
        LinePoint2(1,1,j)=SetOut2(1,1,j);%��ʼ�������
        LinePoint2(1,2,j)=SetOut2(1,NumFlag2(j),j);%�����������
        LinePoint2(2,1,j)=Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j);%��ʼ��������
        LinePoint2(2,2,j)=Result_Ls2(1,j)*SetOut2(1,NumFlag2(j),j)+Result_Ls2(2,j);%������������
    else
%         plot([(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j) ...
%         (SetOut2(2,NumFlag2(j),j)-Result_Ls2(2,j))/Result_Ls2(1,j)],...
%         [SetOut2(2,1,j) SetOut2(2,NumFlag2(j),j)],'k-x','LineWidth',1,'MarkerSize',3);
        LinePoint2(1,1,j)=(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j);%��ʼ�������
        LinePoint2(1,2,j)=(SetOut2(2,NumFlag2(j),j)-Result_Ls2(2,j))/Result_Ls2(1,j);%�����������
        LinePoint2(2,1,j)=SetOut2(2,1,j);%��ʼ��������
%         LinePoint2(2,2,j)=SetOut2(2,NumFlag1(j),j);%������������
        LinePoint2(2,2,j)=SetOut2(2,NumFlag2(j),j);%������������
    end 
end
% hold off

%ʹ��ֱ�ߵĶ˵���ʾֱ��
% figure(104);
% hold on
% for i=1:size(LinePoint1,3)
%    plot([LinePoint1(1,1,i) LinePoint1(1,2,i)],[LinePoint1(2,1,i) LinePoint1(2,2,i)],'g-','LineWidth',2);%����ֱ��
%    %���ƶ˵��ԭ�������
%    plot([LinePoint1(1,1,i) 0],[LinePoint1(2,1,i) 0],'b--');%����ֱ��  
%    plot([0 LinePoint1(1,2,i)],[0 LinePoint1(2,2,i)],'y--');%����ֱ��   
% end
% for i=1:size(LinePoint2,3)
%    plot([LinePoint2(1,1,i) LinePoint2(1,2,i)],[LinePoint2(2,1,i) LinePoint2(2,2,i)],'r-','LineWidth',2);%����ֱ��
%    %���ƶ˵��ԭ�������
%    plot([LinePoint2(1,1,i) 0],[LinePoint2(2,1,i) 0],'r--');%����ֱ��  
%    plot([0 LinePoint2(1,2,i)],[0 LinePoint2(2,2,i)],'m--');%����ֱ��
% end
% hold off


[thetaD1]=LineDomain(LinePoint1,Result_Ls1);
%����ͬ֡���ݵĲ�ֱͬ�߶����򻥳�
if size(thetaD1,2)>2
    for i=1:size(thetaD1,2)-1
%         if (thetaD1(1,i)<=thetaD1(1,i+1)<=thetaD1(2,i))||(thetaD1(1,i+1)<=thetaD1(1,i)<=thetaD1(2,i+1))
        if (thetaD1(1,i)<thetaD1(1,i+1))&&(thetaD1(1,i+1)<thetaD1(2,i))&&(thetaD1(2,i)<thetaD1(2,i+1))
            thetaD1(2,i)=thetaD1(1,i+1)-0.001;%0.001/pi*180=0.057��
        end
    end
end
[thetaD2]=LineDomain(LinePoint2,Result_Ls2);
%����ͬ֡���ݵĲ�ֱͬ�߶����򻥳�
if size(thetaD2,2)>2
    for i=1:size(thetaD2,2)-1
%         if (thetaD1(1,i)<=thetaD1(1,i+1)<=thetaD1(2,i))||(thetaD1(1,i+1)<=thetaD1(1,i)<=thetaD1(2,i+1))
        if (thetaD2(1,i)<thetaD2(1,i+1))&&(thetaD2(1,i+1)<thetaD2(2,i))&&(thetaD2(2,i)<thetaD2(2,i+1))
            thetaD2(2,i)=thetaD2(1,i+1)-0.001;%0.001/pi*180=0.057��
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
            reDomain(1,num)=max(thetaA1,thetaA2);%�غ�ֱ�ߵ���ʼ�ǶȦ�ai
            reDomain(2,num)=min(thetaB1,thetaB2);%�غ�ֱ�ߵĽ����ǶȦ�bi
            reDomain(3,num)=thetaD1(3,i);        %ֱ��1��������x��������ļнǦ�i
            reDomain(4,num)=thetaD1(4,i);        %ԭ�㵽ֱ��1�ľ���li
            reDomain(5,num)=thetaD2(3,j);        %ֱ��2��������x��������ļнǦ�i'
            reDomain(6,num)=thetaD2(4,j);        %ԭ�㵽ֱ��2�ľ���li'
            reDomain(7,num)=thetaD1(5,i);        %ֱ��1�Ĳ���a,(y=ax+b)
            reDomain(8,num)=thetaD1(6,i);        %ֱ��1�Ĳ���b,(y=ax+b)
            reDomain(9,num)=thetaD2(5,j);        %ֱ��2�Ĳ���a,(y=ax+b)
            reDomain(10,num)=thetaD2(6,j);       %ֱ��2,�Ĳ���b,(y=ax+b)
        end
    end
    
end
% figure(105)
% hold on
% for i=1:size(LinePoint1,3)
%    plot([LinePoint1(1,1,i) LinePoint1(1,2,i)],[LinePoint1(2,1,i) LinePoint1(2,2,i)],'g-','LineWidth',2);%����ֱ��    
% end
% for i=1:size(LinePoint2,3)
%    plot([LinePoint2(1,1,i) LinePoint2(1,2,i)],[LinePoint2(2,1,i) LinePoint2(2,2,i)],'r-','LineWidth',2);%����ֱ��
% end
%     for i=1:size(reDomain,2)
%         k1=tan(reDomain(1,i));%����ԭ�㣬����ֱ�߶���������ֱ��б��
%         k2=tan(reDomain(2,i));%����ԭ�㣬����ֱ�߶�����������ֱ��б��
%         %���ֱ��1����������ֱ�ߵĽ���
%         %l1:y=ax+b,k1:y=k1x,k2:y=k2x
%         x11=reDomain(8,i)/(k1-reDomain(7,i));
%         y11=x11*k1;
%         x12=reDomain(8,i)/(k2-reDomain(7,i));
%         y12=x12*k2;
%         
%         %���ֱ��2����������ֱ�ߵĽ���
%         %l1:y=ax+b,k1:y=k1x,k2:y=k2x
%         x21=reDomain(10,i)/(k1-reDomain(9,i));
%         y21=x21*k1;
%         x22=reDomain(10,i)/(k2-reDomain(9,i));
%         y22=x22*k2;
%         
%         plot([0 x11],[0 y11],'m--');%����ֱ��
%         plot([0 x12],[0 y12],'k--');%����ֱ��
%         plot([0 x21],[0 y21],'m--');%����ֱ��
%         plot([0 x22],[0 y22],'k--');%����ֱ��
%         
%         plot([x11 x21],[y11 y21],'m-','LineWidth',2);%����ֱ��
%         plot([x12 x22],[y12 y22],'k-','LineWidth',2);%����ֱ��
%     end
% hold off


sumError=0;

%�ܵ�������
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
%���ϵ��
factor=2*pi/(factor^2);
sumError=sumError*factor;

%����
dx=0;
dy=0;
dTheta=0;
for i=1:size(reDomain,2)
    if reDomain(3,i)==reDomain(5,i)
        %�����dE/d��ʱ��reDomain(3,i)==reDomain(5,i)�ᵼ��dE/d����������󣬵����޷����µ���������
        %��Ϊ��һ��ƫ�ʹ�����߲������
        reDomain(5,i)=reDomain(5,i)+1/180*pi;
        %�����x�������꣩�ĵ���
        dxAi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(1,i)-reDomain(3,i)));
        dxBi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(2,i)-reDomain(3,i)));
        dxTemp=dxBi-dxAi;
        dx=dx+dxTemp;
        %�����y�������꣩�ĵ���
        dyAi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(1,i)-reDomain(3,i)));
        dyBi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             tan(reDomain(2,i)-reDomain(3,i)));
        dyTemp=dyBi-dyAi;
        dy=dy+dyTemp;
        %�����theta����ת�Ƕȣ��ĵ���
        dThetaAi=-reDomain(6,i)^2/(cos(reDomain(1,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(1,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*tan(reDomain(1,i)-reDomain(3,i)));
        dThetaBi=-reDomain(6,i)^2/(cos(reDomain(2,i)-reDomain(5,i))^2)+2*reDomain(4,i)*reDomain(6,i)*(...
                 tan(reDomain(2,i)-reDomain(5,i))/sin(reDomain(3,i)-reDomain(5,i))-...
                 cot(reDomain(3,i)-reDomain(5,i))*tan(reDomain(2,i)-reDomain(3,i)));
        dThetaTemp=dThetaBi-dThetaAi;
        dTheta=dTheta+dThetaTemp;
    else
        %�����x�������꣩�ĵ���
        dxAi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dxBi=2*cos(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dxTemp=dxBi-dxAi;
        dx=dx+dxTemp;
        %�����y�������꣩�ĵ���
        dyAi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(1,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(1,i)-reDomain(3,i))/cos(reDomain(1,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dyBi=2*sin(reDomain(5,i))*(reDomain(6,i)*tan(reDomain(2,i)-reDomain(5,i))-reDomain(4,i)*...
             log(cos(reDomain(2,i)-reDomain(3,i))/cos(reDomain(2,i)-reDomain(5,i)))/...
             sin(reDomain(3,i)-reDomain(5,i)));
        dyTemp=dyBi-dyAi;
        dy=dy+dyTemp;
        %�����theta����ת�Ƕȣ��ĵ���
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