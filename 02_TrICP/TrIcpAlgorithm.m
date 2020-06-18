
%*********************************************************************
%*********************************************************************
%*********************************************************************
%���ܣ� TrICP�㷨�Ķ���ʵ�֣������ֳ���ʵ�ּ���M�ļ�
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
%Step1������ǰ�����ݻ�ȡ�Ͳ�������
%*********************************************************************
%Step1.1�����ü���ɨ���ǵ�ǰ������λ�ã��Ͳ�������
xPre=310;           %����ɨ����ǰһʱ��λ�ú�����
yPre=300;           %����ɨ����ǰһʱ��λ��������
directionPre=0;     %����ɨ����ǰһʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
xLate=330;          %����ɨ���ǵ�ǰʱ��λ�ú�����
yLate=280;          %����ɨ���ǵ�ǰʱ��λ��������
directionLate=20;    %����ɨ���ǵ�ǰʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
numSample=1000;     %���Ƽ���ɨ���ǲ���������ɨ��һȦ�Ĳ�������Ϊ samplingRange/180*samples
samplingRange=140;  %���Ƽ���ɨ���ǲ�����Χ��������ΧΪ(-samplingRange,samplingRange)
MapName='Map13.bmp';%ѡ��ʵ���ͼ��������
ExpData1FromCsv='Exp4_X_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
ExpData2FromCsv='Exp4_P_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
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
title('��Ӹ�˹������ļ������ݵ㼯');    %����
xlabel('ˮƽ����(cm)');     %x��
ylabel('��ֱ����(cm)');     %y��
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'rx','MarkerSize',4); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'gx','MarkerSize',4);  
end
hold off


%*********************************************************************
%*********************************************************************
%Step2��ICP�����㷨
%*********************************************************************
%Step2.1����Ӧ��������
iteraNumber=100;            %�����㷨��������
errorThreshold=0.01;        %������ǰ��ֹ����ѭ���������ֵ
Pii=P;                      %������i��ʱ��P�㼯����ʼֵ��ΪԭʼP
Xi=X;
index=0;
FinalX=0;                   %��¼������λ��
FinalY=0;                   %��¼������λ��
FinalTheta=0;               %��¼��תλ����
setOut1=4:length(P(1,:));
distanceIndiv=zeros(size(P(1,:)));
%*********************************************************************
%Step2.2������ѭ��

for i=1:iteraNumber
    %Step2.2.1����X�㼯��Ѱ��Pi�㼯��ÿ����������,���õ���Ӧ���
    setIn=[Pii;X];                       %����2��n����(Pi��X)����һ��4��n�������(setIn)��������ClosetPointMatch�������Ҫ��
    [setOut0]=ClosetPointMatch1(setIn);  %��ȡ�㼯Pi������㼯
    P0=[setOut0(1,:);setOut0(2,:)];
    X0=[setOut0(3,:);setOut0(4,:)];
    %Step2.2.2��%TrICP���ԣ�������������ŷʽ������С��ǰNp%��������
    [setOut1]=TrStrategy(setOut0,0.7);
    
    Pi=[setOut1(1,:);setOut1(2,:)];
    Xi=[setOut1(3,:);setOut1(4,:)];
    %Step2.2.3������X��Pi����� 
    %coordinateOffste=Xi-Pi; %����Xi��Pi�ĺ��������ֵ
    coordinateOffste=X0-P0; %ԭʼ�����Ե����
    for j=1:length(coordinateOffste(1,:))
        distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%����ÿ��P���X��֮���ŷʽ����
    end%for j=1:length(coordinateOffste(1,:)) 
    %��¼����ת�ǶȺ�λ�����ڵ����еı仯�������ӻ���������ʹ�ã�
    index=index+1;
    error(index)=sum(distanceIndiv)/length(distanceIndiv);%P����X����ƽ�����
    xMatrix(index)=FinalX;              %������λ��
    yMatrix(index)=FinalY;              %������λ��
    thetaMatrix(index)=FinalTheta;      %��ת�Ƕ�
    %Step2.2.4���ж��Ƿ���ǰ��ֹѭ�������error<errorThreshold,����ֹѭ��
    if error(index)<errorThreshold
        break;
    end
    %Step2.2.5�����ݵõ��Ķ�Ӧ��Լ���任����R��T��theta
    [R,T,theta]=CalculateTranMatrix(setOut1);
    FinalX=FinalX+T(1,1);                   %������λ��
    FinalY=FinalY+T(2,1);                   %������λ��
    FinalTheta=FinalTheta+theta/pi*180;     %��ת�Ƕ�
    %Step2.2.6�����ݵõ���R��T�任��P(i+1)
    for k=1:length(Pii(1,:))
        Pii(:,k)=R*Pii(:,k)+T;
    end
end%for i=1:iteraNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.3�����MSE,����һ�����ʵ��Ϊ��ŷ�Ͼ�������㡣

%Step2.3.1����X�㼯��Ѱ��Pi�㼯��ÿ����������,���õ���Ӧ���
PP=Pii;
XX=X;
setIn=[PP;XX];                       %����2��n����(Pi��X)����һ��4��n�������(setIn)��������ClosetPointMatch�������Ҫ��
[setOut1]=ClosetPointMatch1(setIn);  %��ȡ�㼯Pi������㼯
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.3.2������X��Pi����� 
coordinateOffste=Xi-PP;             %����Xi��Pi�ĺ��������ֵ
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%����ÿ��P���X��֮���ŷʽ����
end%for j=1:length(coordinateOffste(1,:)) 
MSE=sum(distanceIndiv)/length(distanceIndiv);%P����X����ƽ�����

%*********************************************************************
%*********************************************************************
%Step3�����ݿ��ӻ�

figure(103);
title('�������');%����
xlabel('��������');%x��
ylabel('�������仯ֵ');%y��
hold on
for i=1:size(error,2)-1
    plot([i i+1],[error(i) error(i+1)],'r-x','MarkerSize',3,'LineWidth',3);  
end
hold off


figure(104);
title('��������仯ֵ');%����
xlabel('��������');%x��
ylabel('��/����������仯ֵ');%y��
hold on
for i=1:size(xMatrix,2)-1
    plot([i i+1],[xMatrix(i) xMatrix(i+1)],'g-x','MarkerSize',3,'LineWidth',3);  
    plot([i i+1],[yMatrix(i) yMatrix(i+1)],'m-x','MarkerSize',3,'LineWidth',3);
end
plot([1 size(xMatrix,2)],[xLate-xPre xLate-xPre],'b-','MarkerSize',3,'LineWidth',3);%������������
plot([1 size(yMatrix,2)],[yPre-yLate yPre-yLate],'r-','MarkerSize',3,'LineWidth',3);%������������
hold off

figure(105)
title('��ת�Ƕȵ����仯ֵ');%����
xlabel('��������');%x��
ylabel('��ת�Ƕȵ����仯ֵ');%y�� 
hold on
for i=1:size(thetaMatrix,2)-1
    plot([i i+1],[thetaMatrix(i) thetaMatrix(i+1)],'k-x','MarkerSize',3,'LineWidth',3);
end
plot([1 size(thetaMatrix,2)],[-directionLate -directionLate],'r-','LineWidth',3);%��ת�Ƕ�������

hold off



figure(106)
title('ƥ��Ч��ͼ');%����
xlabel('ˮƽ����(cm)');%x��
ylabel('��ֱ����(cm)');%y��
% set(gcf,'color','none')
% set(gca,'xtick',[],'ytick',[]) % ͬʱȥ��x��y��Ŀ̶�
hold on
for i=1:size(Pii,2)
    plot([Pii(1,i) Pii(1,i)],[Pii(2,i) Pii(2,i)],'r.','MarkerSize',4);
end
for i=1:size(X,2)
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);
end
hold off

FinalTheta=-FinalTheta;
