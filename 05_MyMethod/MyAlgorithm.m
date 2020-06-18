%*********************************************************************
%*********************************************************************
%*********************************************************************
%���ܣ������򣬸����ֳ���ʵ�ּ���M�ļ�
%���ߣ�Shaofeng Wu 
%ʱ�䣺2019.12.07
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
clear
clc
%*********************************************************************
%*********************************************************************
%Step1�����ݻ�ȡ������Ԥ����
%*********************************************************************
%Step1.1�����ü���ɨ���ǵ�ǰ������λ�ã��Ͳ�������
xPre=310;           %����ɨ����ǰһʱ��λ�ú�����
yPre=300;           %����ɨ����ǰһʱ��λ��������
directionPre=0;     %����ɨ����ǰһʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
xLate=330;          %����ɨ���ǵ�ǰʱ��λ�ú�����
yLate=280;          %����ɨ���ǵ�ǰʱ��λ��������
directionLate=20;   %����ɨ���ǵ�ǰʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
numSample=100;     %���Ƽ���ɨ���ǲ���������ɨ��һȦ�Ĳ�������Ϊ samplingRange/180*samples
samplingRange=180;  %���Ƽ���ɨ���ǲ�����Χ��������ΧΪ(-samplingRange,samplingRange)
MapName='Map32.bmp';%ѡ��ʵ���ͼ��������
ExpData1FromCsv='Exp1_X_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
ExpData2FromCsv='Exp1_P_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
DataSource=false;   %ѡ��������Դģʽ
%*********************************************************************
% Step1.2����ȡ����ɨ����ǰ������λ�õ������Ӧ��X��P�㼯
if DataSource   %ѡ��������Դ��trueģʽΪ�������ɺ���������������
                %ѡ��������Դ��falseģʽΪ�Ӷ�Ӧ��.csv�ļ��ж�ȡ����
                %��Ϊ��ӵĸ�˹����������ģ�ÿ�����ɵ�������һ����
                %����Ϊ�˱�������ʵ�����ݵ�һ���ԣ���Ҫ�����ݱ���������
                %�Է������Ա�ʵ�顣
[setOutX0]=GetLaserDataPointSet(xPre,yPre,numSample,directionPre,11,MapName,samplingRange);   %��ʼλ����Ϊ��xPre,yPre�������ص����ݵ㼯��Ϊģ��㼯X
[setOutP0]=GetLaserDataPointSet(xLate,yLate,numSample,directionLate,22,MapName,samplingRange);%�ƶ�λ����Ϊ��xLate,yLate�������ص����ݵ㼯��Ϊ��ƥ��㼯P

X0=CoordinateTran(setOutX0);%����任��������任Ϊ�ѿ�������
P0=CoordinateTran(setOutP0);%����任��������任Ϊ�ѿ�������
 %���ӻ��������ǲɼ���������(�Ӹ�˹����ǰ)
figure(99);
title('ԭʼ�������ݵ㼯');    %����
hold on
for i=1:size(P0,2)
    plot([P0(1,i) P0(1,i)],[P0(2,i) P0(2,i)],'ro','MarkerSize',4); 
    plot([X0(1,i) X0(1,i)],[X0(2,i) X0(2,i)],'g.','MarkerSize',4);  
end
hold off
%*********************************************************************
% Step1.3�����������Ϊ5dB�ĸ�˹������
[setOutX]=NoiseAddition(setOutX0,5);
[setOutP]=NoiseAddition(setOutP0,5);

%���������ݴ�Ϊ.csv�ļ�
%dlmwrite���Կ��Ʊ���ΪCSV�ļ������ݾ���λ��
dlmwrite('Exp12_X_NoiseAdd_5dB11.csv', setOutX,'precision',32);
dlmwrite('Exp12_P_NoiseAdd_5dB11.csv', setOutP,'precision',32);
else
setOutX = csvread(ExpData1FromCsv);
setOutP = csvread(ExpData2FromCsv);
end %if true

%*********************************************************************
%Step1.4������ת�����������ɼ�����ϵ�任Ϊ�ѿ�������ϵ
X=CoordinateTran(setOutX);%����任��������任Ϊ�ѿ�������
P=CoordinateTran(setOutP);%����任��������任Ϊ�ѿ�������

%���ӻ��������ǲɼ���������(���˸�˹������)
figure(101);
title('��Ӹ�˹������ļ������ݵ㼯');   %����
xlabel('ˮƽ����(cm)');                  %x��
ylabel('��ֱ����(cm)');                  %y��
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off
%*********************************************************************
%*********************************************************************
%Step2�������������ݵ�ֱ�߶���ȡ
%*********************************************************************

%*********************************************************************
%Step2.1�����ݼ�����Ӧ�ָ�
% [SetOut11]=SplitData(X,8,5);%���ݵ�ָ�
% [SetOut22]=SplitData(P,8,5);%���ݵ�ָ�
[SetOut11]=SplitData(X,2,5);%���ݵ�ָ�
[SetOut22]=SplitData(P,2,5);%���ݵ�ָ�
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
    plot([0 0],[0 0],'r.','Markersize',20) %����ԭ�㣨��������λ�ã�
hold off


%*********************************************************************
%Step2.2���Էָ��������Ӽ����ж����Ƿ�Ϊֱ�߶�����������������ٴηָ�

% [SetOut1]=SplitData1(SetOut11,15,3);%�Է�ֱ���������ݼ�ֻ�ָ�һ��
% [SetOut2]=SplitData1(SetOut22,15,3);%�Է�ֱ���������ݼ�ֻ�ָ�һ��
[SetOut1]=SplitData2(SetOut11,15,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�
[SetOut2]=SplitData2(SetOut22,15,3);%�Է�ֱ���������ݼ��ָ�n�Σ�ֱ�������Ӽ���Ϊֱ���������ݼ�

%��ʾ���ݼ��ָ�ʹ��������
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
    plot([0 0],[0 0],'r.','Markersize',20)%����ԭ�㣨��������λ�ã�
hold off

%*********************************************************************
%Step2.3��ֱ�߶����
% [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1);%������С���˷����ֱ�ߣ�Ч���ܺ�
% [Result_Ls2,NumFlag2]=LeastSquareLine1(SetOut2);%������С���˷����ֱ�ߣ�Ч���ܺ�
% [Result_Ls1,NumFlag1]=LineFittingRANSAC(SetOut1);%����RANSAC�����ֱ�ߣ�Ч��һ��
% [Result_Ls2,NumFlag2]=LineFittingRANSAC(SetOut2);%����RANSAC�����ֱ�ߣ�Ч��һ��
[Result_Ls1,NumFlag1]=LineFitting(SetOut1);%RANSAC���Ľ�Ϊ�������ݼ������ֱ�ߣ�Ч���ܺ�
[Result_Ls2,NumFlag2]=LineFitting(SetOut2);%RANSAC���Ľ�Ϊ�������ݼ������ֱ�ߣ�Ч���ܺ�

%��ʾ���յ�������
LinePoint=zeros(2,2,size(SetOut1,3));%�ռ�ֱ�߶εĶ˵�
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
        LinePoint(1,1,j)=SetOut1(1,1,j);%��ʼ�������
        LinePoint(1,2,j)=SetOut1(1,NumFlag1(j),j);%�����������
        LinePoint(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%��ʼ��������
        LinePoint(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%������������
    else
        plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
        (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
        [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-v','LineWidth',1,'MarkerSize',3);
        LinePoint(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%��ʼ�������
        LinePoint(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%�����������
        LinePoint(2,1,j)=SetOut1(2,1,j);%��ʼ��������
        LinePoint(2,2,j)=SetOut1(2,NumFlag1(j),j);%������������
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
%Step3������ֱ�߶���ת�Ƕ�
%*********************************************************************

%*********************************************************************
%Step3.1������ֱ�߲���y=ax+b������ֱ����б�Ƕ�
for i=1:size(Result_Ls1,2)
    dist1(i)=abs(Result_Ls1(2,i)/sqrt(Result_Ls1(1,i)^2+1));
    theta1(i)=atan(Result_Ls1(1,i));
end
for i=1:size(Result_Ls2,2)
    dist2(i)=abs(Result_Ls2(2,i)/sqrt(Result_Ls2(1,i)^2+1));
    theta2(i)=atan(Result_Ls2(1,i));
end    
%*********************************************************************
%Step3.2���Ƕ����˳�㽫�Ƕ��ɻ�����תΪ�Ƕ���
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
%Step3.3��ʹ���ܶȾ���������ǶȲ�
figure(105);
hold on
[DbscanNum,DbscanCell] = DBSCAN(DiffTheta,1,2);%�ܶȾ���
hold off
%ѡȡԪ�����Ĵ�Ϊ�ؼ���
SizeCluster=zeros(1,2);
tempCell=DbscanCell{1,1};
SizeCluster(1,1)=size(tempCell,2);%��¼��ǰ�ذ��������ݸ���
SizeCluster(1,2)=1;               %��¼��Ӧ�ı��
%�����ؼ���
for i=2:DbscanNum
    tempCell=DbscanCell{i,1};
    if size(tempCell,2)>SizeCluster(1,1)%�����ǰ��Ԫ�ش���ǰһ���������
        SizeCluster(1,1)=size(tempCell,2);%��¼��ǰ�ذ��������ݸ���
        SizeCluster(1,2)=i;%��¼��Ӧ�ı��
    end
end
%�����ؼ���
FlagDBC=DbscanCell{SizeCluster(1,2),1};
FinalTheta=0;
%�Թؼ��ص�����������Ȩƽ�����ɵõ�������ת�Ƕ�
for i=1:size(FlagDBC,2)
    FinalTheta=FinalTheta+DiffTheta(1,FlagDBC(1,i));
end
FinalTheta=FinalTheta/size(FlagDBC,2);
% FinalTheta(2,1)=FinalTheta(2,1)/size(FlagDBC,2);

%*********************************************************************
%Step3.4����P�㼯���е���ݵõ��ĽǶ���ת����,Ϊ����ļ���λ����׼��
%�õ���ת�Ƕ�
TranMatrix=[cos(-FinalTheta/180*pi) -sin(-FinalTheta/180*pi);...
            sin(-FinalTheta/180*pi) cos(-FinalTheta/180*pi)];
SetOut2Tra=SetOut2;
%��P�㼯���е���ת
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
%���Ƴ���ת֮���ֱ��
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

%������ѧ��ϵ������תǰ��ֱ�߲����������ת���ֱ�߲�������ʵҲ���Զ���ת��ĵ㼯
%���ֱ�ߣ�����Ϊ�˱���ֱ�߲�����һ���ԣ����Բ�����ѧ��ϵת������
Line2StartPoint=zeros(2,size(SetOut2,3));
for j=1:size(SetOut2,3)
    if abs(SetOut2(1,1,j)-SetOut2(1,NumFlag2(j),j))>5
        %ÿ��ֱ��ѡȡ��ʼ�㣬�������������ֱ�߷������
        Line2StartPoint(1,j)=SetOut2(1,1,j);%��ʼ�������
        Line2StartPoint(2,j)=Result_Ls2(1,j)*SetOut2(1,1,j)+Result_Ls2(2,j);%��ʼ��������
    else
        %ÿ��ֱ��ѡȡ��ʼ�㣬�������������ֱ�߷������
        Line2StartPoint(1,j)=(SetOut2(2,1,j)-Result_Ls2(2,j))/Result_Ls2(1,j);%��ʼ�������
        Line2StartPoint(2,j)=SetOut2(2,1,j);%��ʼ��������
    end 
end

%�Ѿ�֪����ת��ֱ�ߵ�б��Result_Ls2Tra����ת����TranMatrix1
%�Լ�ֱ���ϵ�һ����ʼ��Line2StartPoint
%���Կ��������תǰֱ�ߵ�б��K�ͽؾ�
Result_Ls2Tra=zeros(2,size(Result_Ls2,2));
for i=1:size(Result_Ls2Tra,2)
     Result_Ls2Tra(1,i)=tan(theta2(i)+FinalTheta/180*pi);
end
for i=1:size(Result_Ls2,2)
    temp=Line2StartPoint(:,i)'*(TranMatrix);%��ֱ���ϵĵ���ת
    Result_Ls2Tra(2,i)=temp(2)-Result_Ls2Tra(1,i)*temp(1);   
end
figure(104);
hold on
%��ʾ��ת���ֱ��������
Line1Point=zeros(2,2,size(SetOut1,3));%�ռ�ֱ�߶εĶ˵�
Line2TraPoint=zeros(2,2,size(SetOut2Tra,3));%�ռ�ֱ�߶εĶ˵�
for j=1:size(SetOut1,3)
    if abs(SetOut1(1,1,j)-SetOut1(1,NumFlag1(j),j))>5
        plot([SetOut1(1,1,j) SetOut1(1,NumFlag1(j),j)],...
        [Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j) ...
        Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j)],'b-v','LineWidth',1,'MarkerSize',3);
        Line1Point(1,1,j)=SetOut1(1,1,j);%��ʼ�������
        Line1Point(1,2,j)=SetOut1(1,NumFlag1(j),j);%�����������
        Line1Point(2,1,j)=Result_Ls1(1,j)*SetOut1(1,1,j)+Result_Ls1(2,j);%��ʼ��������
        Line1Point(2,2,j)=Result_Ls1(1,j)*SetOut1(1,NumFlag1(j),j)+Result_Ls1(2,j);%������������
    else
        plot([(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j) ...
        (SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j)],...
        [SetOut1(2,1,j) SetOut1(2,NumFlag1(j),j)],'b-v','LineWidth',1,'MarkerSize',3);
        Line1Point(1,1,j)=(SetOut1(2,1,j)-Result_Ls1(2,j))/Result_Ls1(1,j);%��ʼ�������
        Line1Point(1,2,j)=(SetOut1(2,NumFlag1(j),j)-Result_Ls1(2,j))/Result_Ls1(1,j);%�����������
        Line1Point(2,1,j)=SetOut1(2,1,j);%��ʼ��������
        Line1Point(2,2,j)=SetOut1(2,NumFlag1(j),j);%������������
    end 
end
for j=1:size(SetOut2Tra,3)
    if abs(SetOut2Tra(1,1,j)-SetOut2Tra(1,NumFlag2(j),j))>5
        plot([SetOut2Tra(1,1,j) SetOut2Tra(1,NumFlag2(j),j)],...
        [Result_Ls2Tra(1,j)*SetOut2Tra(1,1,j)+Result_Ls2Tra(2,j) ...
        Result_Ls2Tra(1,j)*SetOut2Tra(1,NumFlag2(j),j)+Result_Ls2Tra(2,j)],'k-o','LineWidth',1,'MarkerSize',3);
        Line2TraPoint(1,1,j)=SetOut2Tra(1,1,j);%��ʼ�������
        Line2TraPoint(1,2,j)=SetOut2Tra(1,NumFlag2(j),j);%�����������
        Line2TraPoint(2,1,j)=Result_Ls2Tra(1,j)*SetOut2Tra(1,1,j)+Result_Ls2Tra(2,j);%��ʼ��������
        Line2TraPoint(2,2,j)=Result_Ls2Tra(1,j)*SetOut2Tra(1,NumFlag2(j),j)+Result_Ls2Tra(2,j);%������������
    else
        plot([(SetOut2Tra(2,1,j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j) ...
        (SetOut2Tra(2,NumFlag2(j),j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j)],...
        [SetOut2Tra(2,1,j) SetOut2Tra(2,NumFlag2(j),j)],'k-o','LineWidth',1,'MarkerSize',3);
        Line2TraPoint(1,1,j)=(SetOut2Tra(2,1,j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j);%��ʼ�������
        Line2TraPoint(1,2,j)=(SetOut2Tra(2,NumFlag2(j),j)-Result_Ls2Tra(2,j))/Result_Ls2Tra(1,j);%�����������
        Line2TraPoint(2,1,j)=SetOut2Tra(2,1,j);%��ʼ��������
        Line2TraPoint(2,2,j)=SetOut2Tra(2,NumFlag2(j),j);%������������
    end 
end
hold off
for i=1:size(Result_Ls2Tra,2)
%     dist2Tra(i)=abs(Result_Ls2Tra(2,i)/sqrt(Result_Ls2Tra(1,i)^2+1));%ԭ�㵽ֱ�߾���
    dist2Tra(i)=dist2(i);%��Ϊֱ����תֻ�ı�ȣ����ı�ѣ�����dist2Tra(i)=abs(Result_Ls2Tra(2,i)/sqrt(Result_Ls2Tra(1,i)^2+1))Ҳ��ȷ�����������ϵ�Ч��ʵ�ʽ�����ȣ�
    theta2Tra(i)=atan(Result_Ls2Tra(1,i));
end 


%*********************************************************************
%*********************************************************************
%Step4������ֱ�߶�λ��
%*********************************************************************


%*********************************************************************
%Step4.1������ֱ�߼�ľ��룬��Ϊ�����ϴ󲿷�ֱ�߶��ǲ�ƽ�еģ�������ƽ�У�
%         ����ֱ�߾������ԭ���ǣ�ȡһ��ֱ�߶ε��е㣬Ȼ�������е㵽����ֱ
%         �߶εľ��롣

%����ÿ���߶ε��е�
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

%��ѡֱ�߶ζ�
Add=1;
CorrectLinePairsFlag=0;
for i=1:size(theta1,2) 
    for j=1:size(theta2Tra,2)
        if abs(theta1(i)-theta2Tra(j))<2/180*pi%2��
            CorrectLinePairsFlag(1,Add)=i;
            CorrectLinePairsFlag(2,Add)=j;
            Add=Add+1;
        end
    end
end




%����ֱ�߶Լ������
DistLines=zeros(1,size(CorrectLinePairsFlag,2));
for i=1:size(CorrectLinePairsFlag,2)
    xMid=Line1CenterPoints(1,CorrectLinePairsFlag(1,i));
    yMid=Line1CenterPoints(2,CorrectLinePairsFlag(1,i));
    DistLines(i)=abs((Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))*xMid-yMid+Result_Ls2Tra(2,CorrectLinePairsFlag(2,i)))/sqrt(Result_Ls2Tra(1,CorrectLinePairsFlag(2,i))^2+1));
end


%*********************************************************************
%Step4.2���Ƕ�У���󣬵㼯P��ֱ�߶εķ�����


%�Ƕ�У���󣬵㼯P��ֱ�߶εķ�����
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
%Step4.3�������е�ֱ�߶Էָ������ֱ�߶�Ϊһ��Ĵ󼯺�
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
        if abs(theta1(LineMatrix(1,1,i))-theta1(LineMatrix(1,1,j)))<5/180*pi%����ֱ�ߵļнǴ����趨��5�ȣ�����õ��������߶ζ�Ϊƽ��ֱ��
            continue;%����н�С����ֵ��������ǰ����
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
%Step4.4����С���˷����λ��

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
%Step4.5��ʹ���ܶȾ������Ż�λ��
figure(106);
hold on
[DbscanNum1,DbscanCell1] = DBSCAN(FinalVectorMatrix,1.5,2);
hold off
%ѡȡԪ�����Ĵ�Ϊ�ؼ���
SizeCluster1=zeros(1,2);
tempCell1=DbscanCell1{1,1};
SizeCluster1(1,1)=size(tempCell1,2);%��¼��ǰ�ذ��������ݸ���
SizeCluster1(1,2)=1;                %��¼��Ӧ�ı��
%�����ؼ���
for i=2:DbscanNum1
    tempCell1=DbscanCell1{i,1};
    if size(tempCell1,2)>SizeCluster1(1,1)%�����ǰ��Ԫ�ش���ǰһ���������
        SizeCluster1(1,1)=size(tempCell1,2);%��¼��ǰ�ذ��������ݸ���
        SizeCluster1(1,2)=i;%��¼��Ӧ�ı��
    end
end
%�����ؼ���
FlagDBC1=DbscanCell1{SizeCluster1(1,2),1};
FinalVector=zeros(2,1);
%�Թؼ��ص�����������Ȩƽ�����ɵõ�������ת�Ƕ�
for i=1:size(FlagDBC1,2)
    FinalVector(1,1)=FinalVector(1,1)+FinalVectorMatrix(1,FlagDBC1(1,i));
    FinalVector(2,1)=FinalVector(2,1)+FinalVectorMatrix(2,FlagDBC1(1,i));
end
FinalVector(1,1)=FinalVector(1,1)/size(FlagDBC1,2);
FinalVector(2,1)=FinalVector(2,1)/size(FlagDBC1,2);


figure(107)
title('ƥ��Ч��ͼ');%����
xlabel('ˮƽ����(cm)');%x��
ylabel('��ֱ����(cm)');%y��
%set(gca,'xtick',[],'ytick',[]) % ͬʱȥ��x��y��Ŀ̶�
hold on
TranMatrix1=[cos(FinalTheta/180*pi) -sin(FinalTheta/180*pi);...
            sin(FinalTheta/180*pi) cos(FinalTheta/180*pi)];
for k=1:length(P(1,:))%����ת��
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
%Step4.6������MSE,����һ�����ʵ��Ϊ��ŷ�Ͼ�������㡣

%Step4.6.1����X�㼯��Ѱ��Pi�㼯��ÿ����������,���õ���Ӧ���
PP=P;
XX=X;
setIn=[PP;XX];                      %����2��n����(Pi��X)����һ��4��n�������(setIn)��������ClosetPointMatch�������Ҫ��
[setOut1]=ClosetPointMatch1(setIn);  %��ȡ�㼯Pi������㼯
Xi=[setOut1(3,:);setOut1(4,:)];
%Step4.6.2������XX��PP����� 
coordinateOffste=Xi-PP;             %����Xi��Pi�ĺ��������ֵ
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%����ÿ��P���X��֮���ŷʽ����
end%for j=1:length(coordinateOffste(1,:)) 

MSE=sum(distanceIndiv)/length(distanceIndiv);%P����X����ƽ�����
    

