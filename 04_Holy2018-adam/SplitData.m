function [SetOut]=SplitData(SetIn,NumPoints,Factor)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：输入激光数据点集，然后按照相邻点间的欧氏距离自动分割成不同的子集
%输入：1）SetIn:待分割的激光数据点集
%         2×length（SetIn）矩阵SetIn=|x|横坐标
%                                     |y|纵坐标
%      2）1×l的Threshold：欧氏距离阈值
%      3）1×lNumPoints:分割得到的子集中数据点数量阈值，数量小于该阈值时，该
%         点集抛弃                     
%输出：1）SetOut:分割后得到的子集的集合，为2×m×n的三维矩阵，n代表2×m的二
%         维矩阵个数，2为二维矩阵的行数，m为二维矩阵的列数。
%         输出的第i个2×m的二维矩阵表示为：
%         SetOut[:,:,i]=|x1 x2 ... xm| 
%                       |y1 y2 ... ym|
%作者：Shaofeng Wu 
%时间：2018.08.09
%邮箱：shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
SetLength=size(SetIn,2);
NumSubSet=1;
FlagList=zeros(1,SetLength);
Threshold=0;
for i=1:SetLength-1
    eDist=sqrt((SetIn(1,i)-SetIn(1,i+1))^2+(SetIn(2,i)-SetIn(2,i+1))^2);
    polarR=sqrt(SetIn(1,i)^2+SetIn(2,i)^2);%极半径
    Threshold=Factor*polarR*sqrt(2-2*cos( 360/681/180*pi));
    if eDist>Threshold
        if NumSubSet>NumPoints%子集点数大于阈值
            FlagList(1,i+1)=1;
        else                  %子集点数小于阈值
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
        SetOut(:,:,SubSetIndex)=0;%全部置零
        j=1;
    end
    SetOut(:,j,SubSetIndex)=SetIn(:,i);
    j=j+1;
end



