

function [R,T,theta]=CalculateTranMatrix1(set,p)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ����������㼯P,X֮�����ת����R����ƽ�ƾ���T��
%���룺��������һһ��Ӧ�ĵ㼯�Եľ���set���Ͳ���p
%       1.setΪ4��n����set(:,i)=| xP  |    
%                                 | yP  |     
%                                 | xX  |    
%                                 | yX  | 
%       2.pΪ����
%���:�������� ��ת����R����ƽ�ƾ���T��
%���У� 1.RΪ2��2����R=|cos(alpha)   sin(alpha)|
%                        |-sin(alpha)  cos(alpha)|
%       2.TΪ2��1����T=|��x|
%                        |��y|
%����(xP,yP)Ϊ��ƥ��㼯P��(xX,yX)ģ��㼯X
%set�㼯���ݻ��ڵѿ�������ϵ
%������λ��Ϊ���ƣ�m��
%���ߣ�Shaofeng Wu 
%ʱ�䣺2017.11.22
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%Step1:��������㼯��Ȩ��
P(1,:)=set(1,:);
P(2,:)=set(2,:);

X(1,:)=set(3,:);
X(2,:)=set(4,:);
Diff=P-X;
for i=1:size(Diff,2)%����ÿ�����ݵăȻ��Ŀ���
	dotDiff(i)=sqrt(dot(Diff(:,i),Diff(:,i)));
end

for j=1:size(dotDiff,2)
    W(j)=p/(dotDiff(j)^(2-p)+1e-8);
end

W_normalized = W/sum(W);


%*********************************************************************
%*********************************************************************
%Step1:��������㼯������
pSetCenterPoint=[0,0];%��ƥ��㼯���ģ�xPCenter,yPCenter��
xSetCenterPoint=[0,0];%ģ��㼯����(xXCenterX,yXCenter)
xP=set(1,:);          %��4��1ά�����ж�ȡ�����������ڼ�������
yP=set(2,:);
xX=set(3,:);
yX=set(4,:);
pSetCenterPoint(1)=xP* W_normalized';  %��ƥ��㼯P���ĺ�����
pSetCenterPoint(2)=yP* W_normalized';   %��ƥ��㼯P����������
xSetCenterPoint(1)=xX* W_normalized';   %ģ��㼯X���ĺ�����
xSetCenterPoint(2)=yX* W_normalized';   %ģ��㼯X����������
% pSetCenterPoint(1)=sum(xP)/numel(xP);   %��ƥ��㼯P���ĺ�����
% pSetCenterPoint(2)=sum(yP)/numel(yP);   %��ƥ��㼯P����������
% xSetCenterPoint(1)=sum(xX)/numel(xX);   %ģ��㼯X���ĺ�����
% xSetCenterPoint(2)=sum(yX)/numel(yX);   %ģ��㼯X����������
%*********************************************************************
%*********************************************************************
%Step2:�ֱ���������㼯P��Xÿ���㵽��㼯���ĵ�ƫ��
xPOffset=xP- pSetCenterPoint(1);
yPOffset=yP- pSetCenterPoint(2);
xXOffset=xX- xSetCenterPoint(1);
yXOffset=yX- xSetCenterPoint(2);
%*********************************************************************
%*********************************************************************
%Step3:����任����R,T
Numerator=sum(xXOffset.*yPOffset-yXOffset.*xPOffset);   %����
Denominator=sum(xXOffset.*xPOffset+yXOffset.*yPOffset); %��ĸ
alpha=atan(Numerator/Denominator);                      %ƫ�ƽǶȽǶ�
R=[cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];       %��ת����
T=(xSetCenterPoint')-R*(pSetCenterPoint');              %ƽ�ƾ��� 
theta=alpha;
 


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 