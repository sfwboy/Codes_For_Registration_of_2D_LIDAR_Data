clc
clear
%*********************************************************************
%Step1.1���Ӵ洢��csv�ļ��ж�ȡ����
setOutX = csvread('Exp4_X_NoiseAdd_5dB.csv');
setOutP = csvread('Exp4_P_NoiseAdd_5dB.csv');
FinalTheta=-26.6738/180*pi;               %����PSMԴ��õ�����ת�Ƕ�
FinalVector=[24.0180,-35.7921]; %����PSMԴ��õ���ƽ�ƾ���
%*********************************************************************
%Step1.2������ת�����������ɼ�����ϵ�任Ϊ�ѿ�������ϵ
X=CoordinateTran(setOutX);%����任��������任Ϊ�ѿ�������
P=CoordinateTran(setOutP);%����任��������任Ϊ�ѿ������� 
%Step1.3����ʾԭʼ����
figure(1111111);
% title('ԭʼ�������ݵ�');    %����
% xlabel('ˮƽ����(cm)');     %x��
% ylabel('��ֱ����(cm)');     %y��
hold on
for i=1:size(P,2)
    plot([P(1,i) P(1,i)],[P(2,i) P(2,i)],'r.','MarkerSize',4); 
     plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off

% FinalVector=[20,-20];
TranMatrix1=[cos(FinalTheta) -sin(FinalTheta);...
            sin(FinalTheta) cos(FinalTheta)];
for k=1:length(P(1,:))%����ת��
    P1(:,k)=TranMatrix1*P(:,k);
end
for k=1:length(P1(1,:))
    P1(1,k)=P1(1,k)+FinalVector(1);
    P1(2,k)=P1(2,k)+FinalVector(2);
end

figure(1111112);
% title('ԭʼ�������ݵ�');    %����
% xlabel('ˮƽ����(cm)');     %x��
% ylabel('��ֱ����(cm)');     %y��
hold on
for i=1:size(P1,2)
    plot([P1(1,i) P1(1,i)],[P1(2,i) P1(2,i)],'r.','MarkerSize',4); 
     plot([X(1,i) X(1,i)],[X(2,i) X(2,i)],'g.','MarkerSize',4);  
end
hold off








% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%���MSE,����һ�����ʵ��Ϊ��ŷ�Ͼ�������㡣
% %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step2.2.1����X�㼯��Ѱ��Pi�㼯��ÿ����������,���õ���Ӧ���
PP=P1;
XX=X;
setIn=[PP;XX];                       %����2��n����(Pi��X)����һ��4��n�������(setIn)��������ClosetPointMatch�������Ҫ��
[setOut1]=ClosetPointMatch1(setIn);  %��ȡ�㼯Pi������㼯
Xi=[setOut1(3,:);setOut1(4,:)];
%Step2.2.2������X��Pi����� 
coordinateOffste=Xi-PP;             %����Xi��Pi�ĺ��������ֵ
for j=1:length(coordinateOffste(1,:))
    distanceIndiv(j)=sqrt(coordinateOffste(1,j)^2+coordinateOffste(2,j)^2);%����ÿ��P���X��֮���ŷʽ����
end%for j=1:length(coordinateOffste(1,:)) 

MSE=sum(distanceIndiv)/length(distanceIndiv);%P����X����ƽ�����