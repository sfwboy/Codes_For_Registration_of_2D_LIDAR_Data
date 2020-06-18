function [Result,NumFlag]=LineFitting(DataIn)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%函数功能：输入点集，拟合点集中包含的直线段，输出直线参数
%          对于LineFittingRANSAC()的改进，由随机选取点变为遍历所有点
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
%时间：2018.10.06
%邮箱：shaofeng693@126.com
%识别数据中的非零项
for k=1:size(DataIn,3)
    PointSet=DataIn(:,:,k);
    NonZeroNum=0;
    for i=1:size(PointSet,2)
        if PointSet(1,i)~=0 || PointSet(2,i)~=0
            NonZeroNum=NonZeroNum+1;%计算非零项的个数
        end
    end
    SumDistPointsToLine=300000;
    tempSumDistPointsToLine=0;
    for i=1:NonZeroNum-1
        for m=i+1:NonZeroNum
            x1=PointSet(1,i);
            y1=PointSet(2,i);
            x2=PointSet(1,m);
            y2=PointSet(2,m);
            %从点集中选取两个随机点拟合直线
            slope=(y2-y1)/(x2-x1);
            intercept=y1-slope*x1;
            %计算所有点到该直线的距离之和
            for j=1:NonZeroNum
                tempSumDistPointsToLine=abs(slope*PointSet(1,j)-PointSet(2,j)+intercept)...
                /sqrt(slope^2+1);
            end
            if tempSumDistPointsToLine<SumDistPointsToLine
                SumDistPointsToLine=tempSumDistPointsToLine;
                NumFlag(k)=NonZeroNum;
                Result(1,k)=slope;
                Result(2,k)=intercept;
            end
        end
    end
end
