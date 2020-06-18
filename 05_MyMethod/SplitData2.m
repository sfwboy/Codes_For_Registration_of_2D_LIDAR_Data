function [SetOut]=SplitData2(SetIn,LenThesshold,LineThreshold)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%函数功能：1.根据阈值Threshold对输入的点集SetIn进行直线段的判断，如果是直线
%           数据（小于阈值）则不处理，如果不是直线（曲线，大于阈值）则对数据
%           集分割，直到分割的子集小于阈值（认为是直线）为止。
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
%时间：2018.10.19
%邮箱：shaofeng693@126.com
%**************************************************************************
%**************************************************************************
%**************************************************************************

Flag=true;
NumSetOut=1;
for i=1:size(SetIn,3)
    data=0;%清零
    add0=1;
    for j=1:size(SetIn(:,:,i),2)
        if SetIn(1,j,i)~=0
            data(1,add0,1)=SetIn(1,j,i);
            data(2,add0,1)=SetIn(2,j,i);
            add0=add0+1;
        end
    end
    Add=0;
    %*********************************************************************
    %*********************************************************************
    %Step1：判断当前点集是否为直线段集合，不是则分割，分割点为离首尾两点拟合
    %       的直线最远的点
    while Flag
        for add1=1:size(data,3)
            %*********************************************************************
            %Step1.1：拟合经过当前数据集首尾两点的直线：y=ax+b
            temp=data;
            data=0;
            Add=0;
            for ii=1:size(temp,3)
            len=0;
            for j=1:size(temp(:,:,ii),2)
                if temp(1,j,ii)~=0
                len=len+1;
                end
            end
            %利用两点式计求直线参数
            lineA=(temp(2,1,ii)-temp(2,len,ii))/(temp(1,1,ii)-temp(1,len,ii));%斜率
            lineB=temp(2,1,ii)-temp(1,1,ii)*lineA;                            %截距
            distMax=0;
            maxFlag=0;
            %*********************************************************************
            %Step1.2：计算当前集合所有点到直线的距离，保存最大值和对应的坐标
            for k=1:len
                %点到直线距离公式，注意：y=ax+b与Ax+By+C=0之间的参数变换
                distTemp=abs(lineA*temp(1,k,ii)-temp(2,k,ii)+lineB)/sqrt(lineA^2+1);
                if distMax<distTemp
                    distMax=distTemp;
                    maxFlag=k;
                end   
            end
            %*********************************************************************
            %Step1.3：分割集合
            if distMax>LineThreshold%大于阈值，则分为两个子集
                add2=1;
                Add=Add+1;
                for m=1:maxFlag
                    data(1,add2,Add)=temp(1,m,ii);
                    data(2,add2,Add)=temp(2,m,ii);
                    add2=add2+1;
                end
                add2=1;
                Add=Add+1;
                for n=maxFlag:1:len
                    data(1,add2,Add)=temp(1,n,ii);
                    data(2,add2,Add)=temp(2,n,ii);
                    add2=add2+1;
                end
            else      %小于阈值则不分割
                for jj=1:len
                    data(1,jj,ii)=temp(1,jj,ii);
                    data(2,jj,ii)=temp(2,jj,ii); 
                end
                Add=ii;
            end
            end
        end
        %循环终止条件，当前集合不能再分割，程序标记为temp和data矩阵中包含的二维矩阵个数不一样
    if size(temp,3)~=size(data,3)
        Flag=true;
    else
        Flag=false;
    end
    end
    %*********************************************************************
    %*********************************************************************
    %Step2：判断集合包含直线段长度，大于阈值LenThesshold才输出
    for t=1:size(data,3)
        numNotZero=0;
        for r=1:size(data,2)
            if data(1,r,t)~=0
                numNotZero=numNotZero+1;
            end
        end
        if sqrt((data(1,1,t)-data(1,numNotZero,t))^2+(data(2,1,t)-data(2,numNotZero,t))^2)>LenThesshold
            for k=1:numNotZero
                SetOut(1,k,NumSetOut)=data(1,k,t);
                SetOut(2,k,NumSetOut)=data(2,k,t);
            end
            NumSetOut=NumSetOut+1;
        end
    end
    Flag=true; 
end




