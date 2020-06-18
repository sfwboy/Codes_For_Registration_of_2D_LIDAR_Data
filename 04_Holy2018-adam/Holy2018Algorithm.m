%*********************************************************************
%*********************************************************************
%*********************************************************************
%���ܣ� �����򣬸����ֳ���ʵ�ּ���M�ļ�
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
xPre=250;           %����ɨ����ǰһʱ��λ�ú�����
yPre=250;           %����ɨ����ǰһʱ��λ��������
directionPre=0;     %����ɨ����ǰһʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
xLate=300;          %����ɨ���ǵ�ǰʱ��λ�ú�����
yLate=300;          %����ɨ���ǵ�ǰʱ��λ��������
directionLate=30;    %����ɨ���ǵ�ǰʱ�̳������׼����ļнǣ���λ���ȣ�����׼������Ϊˮƽ����
numSample=1000;     %���Ƽ���ɨ���ǲ���������ɨ��һȦ�Ĳ�������Ϊ samplingRange/180*samples
samplingRange=140;  %���Ƽ���ɨ���ǲ�����Χ��������ΧΪ(-samplingRange,samplingRange)
MapName='Map13.bmp';%ѡ��ʵ���ͼ��������
ExpData1FromCsv='Exp14_X_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
ExpData2FromCsv='Exp14_P_NoiseAdd_5dB.csv';%��ȡ��csv�ļ�����/���ݴ�Ϊcsv�ļ�������
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
% xlabel('ˮƽ����(cm)');     %x��
% ylabel('��ֱ����(cm)');     %y��
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'rx','MarkerSize',3); 
    plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'gx','MarkerSize',3);  
end
hold off

%*********************************************************************
%*********************************************************************
%Step2��Registration of lines�㷨
%*********************************************************************
%Step2.1����Ӧ��������

%adam������Ӧѧϰ���㷨����������
parameters=[0;0;25];    %zeros(3,1);�ֱ��ʾ������λ�ơ�������λ�ƺ���ת�Ƕ�
learningRate = 0.1;     %��ʼѧϰ��
b1 = 0.9;               %adam�㷨���߽����Ĭ��ֵ
b2 = 0.999;             %adam�㷨���߽����Ĭ��ֵ
e = 0.00000001;         %adam�㷨���߽����Ĭ��ֵ
mt=zeros(3,1);          %һ�׶���
vt=zeros(3,1);          %���׶���


thresholdError=100;     %��ֵ


Xk=X;
Pk=P;
% iter=0;
iter1=0;


for i=1:2500
%*********************************************************************
%Step2.2�����ú�������ֱ���ݶ�
[dx,dy,dTheta,sumError]=Registration_of_lines(Xk,Pk);
%*********************************************************************
%Step2.3����������
%adam�Ż��㷨��������
gradient(1,1)=real(dx);       %������λ���ݶ�,ֻȡʵ��
gradient(2,1)=real(dy);       %������λ���ݶȣ�ֻȡʵ��
gradient(3,1)=real(dTheta);   %��ת�Ƕ��ݶȣ�ֻȡʵ��

mt=b1*mt+(1-b1)*gradient;       %����һ�׶���
vt=b2*vt+(1-b2)*(gradient.^2);  %���¶��׶���
mtt=mt/(1-(b1.^(i+1)));         %�������
vtt=vt/(1-(b2.^(i+1)));         %�������
vtt_sqrt=[sqrt(vtt(1,1));sqrt(vtt(2,1));sqrt(vtt(3,1))];
%�ݶ��½��������������
parameters=parameters-learningRate*mtt./(vtt_sqrt+e);
finalDx=parameters(1,1);
finalDy=parameters(2,1);
finalTheta=parameters(3,1);
%*********************************************************************
%Step2.4���������ò����Ե㼯������Ӧ�任
R=[cos(finalTheta/180*pi) -sin(finalTheta/180*pi);sin(finalTheta/180*pi) cos(finalTheta/180*pi)];       %��ת����
T=[finalDx;finalDy];
for j=1:size(P,2)
        Pk(:,j)=R*P(:,j)+T;
end
%*********************************************************************
%Step2.5����¼����ת�ǶȺ�λ�����ڵ����еı仯�������ӻ���������ʹ�ã�
iter1=iter1+1;
res(1,iter1)=parameters(1,1);
res(2,iter1)=parameters(2,1);
res(3,iter1)=parameters(3,1);
error(iter1)=sumError; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.6:���MSE,����һ�����ʵ��Ϊ��ŷ�Ͼ�������㡣

%Step2.6.1����X�㼯��Ѱ��Pi�㼯��ÿ����������,���õ���Ӧ���
PP=Pk;
XX=Xk;
setIn=[PP;XX];                       %����2��n����(Pi��X)����һ��4��n�������(setIn)��������ClosetPointMatch�������Ҫ��
[setOut1]=ClosetPointMatch1(setIn);  %��ȡ�㼯Pi������㼯
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.6.2������X��Pi����� 
coordinateOffste=Xi-PP;             %����Xi��Pi�ĺ��������ֵ
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%����ÿ��P���X��֮���ŷʽ����
end%for j=1:length(coordinateOffste(1,:)) 
MSE=sum(distanceIndiv)/length(distanceIndiv);%P����X����ƽ�����



%*********************************************************************
%*********************************************************************
%Step3�����ݿ��ӻ�
figure(13);
title('�������');%����
xlabel('��������');%x��
ylabel('�������仯ֵ');%y��
hold on
for i=1:size(error,2)-1
    plot([i i+1],[error(i) error(i+1)],'k*-');
end
hold off

figure(14);
title('��������仯ֵ');%����
xlabel('��������');%x��
ylabel('��/����������仯ֵ');%y��
hold on
for i=1:size(res,2)-1
    plot([i i+1],[res(1,i) res(1,i+1)],'g*-');
    plot([i i+1],[xLate-xPre xLate-xPre],'b-','LineWidth',2);%������
    
    plot([i i+1],[res(2,i) res(2,i+1)],'m*-');
    plot([i i+1],[yPre-yLate yPre-yLate],'r-','LineWidth',2);%������
end
hold off

figure(15)
title('��ת�Ƕȵ����仯ֵ');%����
xlabel('��������');%x��
ylabel('��ת�Ƕȵ����仯ֵ');%y�� 
hold on
for i=1:size(res,2)-1
    plot([i i+1],[res(3,i) res(3,i+1)],'k*-');
    plot([i i+1],[directionLate directionLate],'r-','LineWidth',2);%������
end
hold off


figure(16)
title('ƥ��Ч��ͼ');%����
xlabel('ˮƽ����(cm)');%x��
ylabel('��ֱ����(cm)');%y��
hold on
for i=1:size(Pk,2)
    plot([Pk(1,i) Pk(1,i)],[Pk(2,i) Pk(2,i)],'r.','MarkerSize',4);
end
for i=1:size(Xk,2)
    plot([Xk(1,i) Xk(1,i)],[Xk(2,i) Xk(2,i)],'g.','MarkerSize',4);
end
hold off

































