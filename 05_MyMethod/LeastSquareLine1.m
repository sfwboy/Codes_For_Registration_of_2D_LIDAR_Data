function [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%函数功能：输入点集，利用最小二乘法拟合点集中包含的直线段，输出直线参数
%输入：三维数据点集DataIn
%     其中，1）DataIn为2×p×q的三维矩阵，q代表2×p的二维矩阵个数，2为二维矩
%              阵的行数，p为二维矩阵的列数。
%              输入的第j个2×p的二维矩阵表示为：
%              SetOut[:,:,j]=|x1 x2 ... xj|
%                            |y1 y2 ... yj|
%              x,y分别为笛卡尔坐标系中的横纵坐标
%输出：二维数据集数Result_Ls,NumFlag
%     其中，1）Result_Ls为2×q的二维矩阵，2为二维矩阵的行数，q为二维矩阵的列
%              数(代表直线段的个数)。
%              输出的第i个2×q的二维矩阵表示为：
%              Result_Ls[2,i]=|a|
%                             |b|
%             拟合直线段在笛卡尔坐标系中表示为：y=ax+b，a为斜率，b为截距。
%          2）NumFlag为1×q的二维矩阵，q为二维矩阵的列数(代表直线段的个数)。
%             NumFlag=|Num1 Num2 Num3 ... Numq|
%             Num代表对应直线段的对应集合中非零点的个数，主要用于后续程序的
%             绘图。
%作者：Shaofeng Wu 
%时间：2018.08.31
%邮箱：shaofeng693@126.com
%**************************************************************************
%**************************************************************************
%**************************************************************************

%最小二乘法拟合直线
[Result_Ls1,NumFlag1]=LeastSquareLine(SetOut1);
%因为最小二乘法对于夹角接近90度的直线拟合存在误差，所以进行以下操作
%%判断所求直线的斜率，如果大于arctan2=63.43度则将该分割点集旋转90度，
%然后再使用最小二乘法拟合直线，最后再把直线旋转回来
% TranMatrix1=[cos(pi/2) -sin(pi/2);sin(pi/2i) cos(pi/2)];
Data1=SetOut1;
TranMatrix1=[0 1;-1 0];
TranFlag1=zeros(1,size(Result_Ls1,2));
for i=1:size(Result_Ls1,2)
    if abs(Result_Ls1(1,i))>6
        TranFlag1(i)=1;
        tempData=Data1(:,:,i);
        for j=1:size(tempData,2)
            temp=tempData(:,j)'*TranMatrix1;
            Data1(:,j,i)=temp';
        end
    end
end
%拟合旋转后的直线
[Result_Ls1Tra0,NumFlag1Tra]=LeastSquareLine(Data1);

Line1StartPoint=zeros(2,size(Data1,3));
for j=1:size(Data1,3)
    if abs(Data1(1,1,j)-Data1(1,NumFlag1Tra(j),j))>5
        if TranFlag1(j)~=0
             %每条直线选取起始点，用来进行下面的直线方程求解
            Line1StartPoint(1,j)=Data1(1,1,j);%起始点横坐标
            Line1StartPoint(2,j)=Result_Ls1Tra0(1,j)*Data1(1,1,j)+Result_Ls1Tra0(2,j);%起始点纵坐标
        end
    else
        if TranFlag1(j)~=0
             %每条直线选取起始点，用来进行下面的直线方程求解
            Line1StartPoint(1,j)=(Data1(2,1,j)-Result_Ls1Tra0(2,j))/Result_Ls1Tra0(1,j);%起始点横坐标
            Line1StartPoint(2,j)=Data1(2,1,j);%起始点纵坐标
        end
    end 
end

%已经知道旋转后直线的斜率Result_Ls1Tra，Result_Ls2Tra和旋转矩阵TranMatrix1
%以及直线上的一个起始点Line1StartPoint，Line2StartPoint
%所以可以求解旋转前直线的斜率K和截距

%因为旋转90度，所以斜率K1=-1/Result_Ls1Tra
for i=1:size(TranFlag1,2)
    if TranFlag1(i)~=0
        Result_Ls1(1,i)=1/Result_Ls1Tra0(1,i);
        temp=Line1StartPoint(:,i)'*(-TranMatrix1);%将直线上的点旋转
        Result_Ls1(2,i)=temp(2)-Result_Ls1(1,i)*temp(1);
    end
end