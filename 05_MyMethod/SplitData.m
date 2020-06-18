function [SetOut]=SplitData(SetIn,NumPoints,Factor)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ����뼤�����ݵ㼯��Ȼ�������ڵ���ŷ�Ͼ����Զ��ָ�ɲ�ͬ���Ӽ�
%���룺1��SetIn:���ָ�ļ������ݵ㼯
%         2��length��SetIn������SetIn=|x|������
%                                     |y|������
%      2��1��l��Threshold��ŷ�Ͼ�����ֵ
%      3��1��lNumPoints:�ָ�õ����Ӽ������ݵ�������ֵ������С�ڸ���ֵʱ����
%         �㼯����                     
%�����1��SetOut:�ָ��õ����Ӽ��ļ��ϣ�Ϊ2��m��n����ά����n����2��m�Ķ�
%         ά���������2Ϊ��ά�����������mΪ��ά�����������
%         ����ĵ�i��2��m�Ķ�ά�����ʾΪ��
%         SetOut[:,:,i]=|x1 x2 ... xm| 
%                       |y1 y2 ... ym|
%���ߣ�Shaofeng Wu 
%ʱ�䣺2018.08.09
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
SetLength=size(SetIn,2);
NumSubSet=1;
FlagList=zeros(1,SetLength);
Threshold=0;
for i=1:SetLength-1
    eDist=sqrt((SetIn(1,i)-SetIn(1,i+1))^2+(SetIn(2,i)-SetIn(2,i+1))^2);
    polarR=sqrt(SetIn(1,i)^2+SetIn(2,i)^2);%���뾶
    Threshold=Factor*polarR*sqrt(2-2*cos( 360/681/180*pi));
    if eDist>Threshold
        if NumSubSet>NumPoints%�Ӽ�����������ֵ
            FlagList(1,i+1)=1;
        else                  %�Ӽ�����С����ֵ
            FlagList(i+1)=2;
        end
        NumSubSet=1;
    end
     NumSubSet=NumSubSet+1;
end
SubSetIndex=1;
j=1;
for i=1:SetLength-1  
    if FlagList(i)==1
        SubSetIndex=SubSetIndex+1;
        j=1;
    end
    if FlagList(i)==2
        SubSetIndex=SubSetIndex;
        SetOut(:,:,SubSetIndex)=0;%ȫ������
        j=1;
    end
    SetOut(:,j,SubSetIndex)=SetIn(:,i);
    j=j+1;
end



