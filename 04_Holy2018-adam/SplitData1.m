function [SetOut2]=SplitData1(SetIn,LenThesshold,LineThreshold)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：1.根据阈值Threshold对输入的点集SetIn进行直线段的判断，如果是直线
%           数据（小于阈值）则不处理，如果不是直线（曲线，大于阈值）则对数据
%           集分割，只分割一次。
%          2.分割完后得到一堆子集，子集代表的直线段长度大于LineThreshold才作
%          为SetOut输出。
%输入：点集SetIn，长度阈值LenThesshold，直线阈值LineThreshold
%     其中，1）SetIn为2×m×n的三维矩阵，n代表2×m的二维矩阵个数，2为二维矩
%              阵的行数，m为二维矩阵的列数。
%              输入的第i个2×m的二维矩阵表示为：
%              SetIn[:,:,i]=|x1 x2 ... xm|
%                           |y1 y2 ... ym|
%              x,y分别为笛卡尔坐标系中的横纵坐标
%           2）LenThesshold为1×1的数值，代表直线长度阈值
%           3）LineThreshold为1×1的数值，代表直线段判断的阈值
%输出：数据集SetOut
%     其中，1）SetOut为2×p×q的三维矩阵，q代表2×p的二维矩阵个数，2为二维矩
%              阵的行数，p为二维矩阵的列数。
%              输出的第j个2×p的二维矩阵表示为：
%              SetOut[:,:,j]=|x1 x2 ... xj|
%                            |y1 y2 ... yj|
%              x,y分别为笛卡尔坐标系中的横纵坐标
%作者：Shaofeng Wu 
%时间：2018.09.10
%邮箱：shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
SetOut2=0;
SetOut1(:,:,1)=SetIn(:,:,1);
Flag=0;
for d=1:size(SetIn,3)
    len=0;
    data1=SetIn(:,:,d);
    for i=1:size(data1,2)
        if data1(1,i)~=0
            len=len+1;
        end
    end
    %*********************************************************************
    %*********************************************************************
    %Step1：判断当前点集是否为直线段集合，不是则分割，分割点为离首尾两点拟合
    %       的直线最远的点
    Add0=0;
    %*********************************************************************
    %Step1.1：拟合经过当前数据集首尾两点的直线：y=ax+b
    %利用两点式计求直线参数
    lineA=(data1(2,1)-data1(2,len))/(data1(1,1)-data1(1,len));
    lineB=data1(2,1)-data1(1,1)*lineA;
    distMax=0;
    maxFlag=0;
    %*********************************************************************
    %Step1.2：计算当前集合所有点到直线的距离，保存最大值和对应的坐标
    for i=1:len
        %点到直线距离公式，注意：y=ax+b与Ax+By+C=0之间的参数变换
        distTemp=abs(lineA*data1(1,i)-data1(2,i)+lineB)/sqrt(lineA^2+1);
        if distMax<distTemp
            distMax=distTemp;
            maxFlag=i;
        end   
    end
    tempData=0;
    %*********************************************************************
    %Step1.3：分割集合
    if distMax>LineThreshold%大于阈值，则分为两个子集
        add=1;
        for i=1:maxFlag
            tempData(1,add,1)=data1(1,i);
            tempData(2,add,1)=data1(2,i);
            add=add+1;
        end
        add=1;
        for j=maxFlag:1:len
            tempData(1,add,2)=data1(1,j);
            tempData(2,add,2)=data1(2,j);
            add=add+1;
        end
        Add0=2;
    else %小于阈值，则不分割
        for j=1:len
            tempData(1,j,1)=data1(1,j);
            tempData(2,j,1)=data1(2,j);
        end
        Add0=1;
    end
    if Add0==1
        for i=1:size(tempData(:,:,1),2)
            SetOut1(1,i,Flag+1)=tempData(1,i,1);
            SetOut1(2,i,Flag+1)=tempData(2,i,1);
        end
        Flag=Flag+1;
    else
        for i=1:size(tempData(:,:,1),2)
            SetOut1(1,i,Flag+1)=tempData(1,i,1);
            SetOut1(2,i,Flag+1)=tempData(2,i,1);
        end        
        for i=1:size(tempData(:,:,2),2)
            SetOut1(1,i,Flag+2)=tempData(1,i,2);
            SetOut1(2,i,Flag+2)=tempData(2,i,2);
        end       
            Flag=Flag+2;
    end
end
Add0=1;
for i=1:size(SetOut1,3)
    add=0;
    for j=1:size(SetOut1,2)
        if SetOut1(1,j,i)~=0
            add=add+1;
        end
    end
    %*********************************************************************
    %*********************************************************************
    %Step2：判断集合包含直线段长度，大于阈值LenThesshold才输出
    if sqrt((SetOut1(1,1,i)-SetOut1(1,add,i))^2+(SetOut1(2,1,i)-SetOut1(2,add,i))^2)>LenThesshold
        for k=1:add
            SetOut2(1,k,Add0)=SetOut1(1,k,i);
            SetOut2(2,k,Add0)=SetOut1(2,k,i);
        end
        Add0=Add0+1;
    end
end





