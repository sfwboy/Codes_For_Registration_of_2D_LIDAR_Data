

function [setOut]=CoordinateTran(setIn)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ�����任����������ϵ����ת��Ϊ�ѿ�������ϵ����
%���룺����Num���������ĵ㼯setIn
%���������Num���ѿ�������ĵ㼯setOut
%setIn=  |range|         setOut=| x  |
%        |theta|                | y  |             
%       ����setOut���ݻ��ڼ�����ϵ
%       ����setIn���ݻ��ڵѿ�������ϵ
%���ߣ�Shaofeng Wu 
%ʱ�䣺2017.11.28
%*********************************************************************
%*********************************************************************
%*********************************************************************
%*********************************************************************
%*********************************************************************
range=setIn(1,:);
theta=setIn(2,:);
%*********************************************************************
%*********************************************************************
%Step1�������ݵ��ɼ�����ϵת��Ϊ�ѿ�������ϵ
alpha=theta*pi/180;   %�Ƕ���תΪ������

%*********************************************************************
%*********************************************************************
%Step2������ת������
x=zeros(size(theta)); %����һ����С��Thetaһ��������
y=zeros(size(theta)); %����һ����С��Thetaһ��������

for m = 1:size(theta,2)
    x(m)=range(m)*cos(alpha(m));
    y(m)=range(m)*sin(alpha(m));
end

setOut=[x;y];








