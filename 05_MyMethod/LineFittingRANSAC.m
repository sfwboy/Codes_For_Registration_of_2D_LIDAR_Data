function [Result,NumFlag]=LineFittingRANSAC(DataIn)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%�������ܣ�����㼯������RANSAC�㷨��ϵ㼯�а�����ֱ�߶Σ����ֱ�߲���
%���룺��ά���ݵ㼯DataIn
%     ���У�1��DataInΪ2��p��q����ά����q����2��p�Ķ�ά���������2Ϊ��ά��
%              ���������pΪ��ά�����������
%              ����ĵ�j��2��p�Ķ�ά�����ʾΪ��
%              SetOut[:,:,j]=|x1 x2 ... xj|
%                            |y1 y2 ... yj|
%              x,y�ֱ�Ϊ�ѿ�������ϵ�еĺ�������
%�������ά���ݼ���Result_Ls,NumFlag
%     ���У�1��Result_LsΪ2��q�Ķ�ά����2Ϊ��ά�����������qΪ��ά�������
%              ��(����ֱ�߶εĸ���)��
%              ����ĵ�i��2��q�Ķ�ά�����ʾΪ��
%              Result_Ls[2,i]=|a|
%                             |b|
%             ���ֱ�߶��ڵѿ�������ϵ�б�ʾΪ��y=ax+b��aΪб�ʣ�bΪ�ؾࡣ
%          2��NumFlagΪ1��q�Ķ�ά����qΪ��ά���������(����ֱ�߶εĸ���)��
%             NumFlag=|Num1 Num2 Num3 ... Numq|
%             Num�����Ӧֱ�߶εĶ�Ӧ�����з����ĸ�������Ҫ���ں��������
%             ��ͼ��
%���ߣ�Shaofeng Wu 
%ʱ�䣺2018.08.31
%���䣺shaofeng693@126.com
%**************************************************************************
%**************************************************************************
%**************************************************************************
%ʶ�������еķ�����
for k=1:size(DataIn,3)
    PointSet=DataIn(:,:,k);
    NonZeroNum=0;
    for i=1:size(PointSet,2)
        if PointSet(1,i)~=0 || PointSet(2,i)~=0
            NonZeroNum=NonZeroNum+1;%���������ĸ���
        end
    end
    NonZeroNum=NonZeroNum-1;
    Iter=5000;
    SumDistPointsToLine=300000;
    tempSumDistPointsToLine=0;
    FinalSlope=0;
    FinalIntercept=0;
    for i=1:Iter
        RandValue=UserRand(NonZeroNum);%���������
        x1=PointSet(1,RandValue(1));
        y1=PointSet(2,RandValue(1));
        x2=PointSet(1,RandValue(2));
        y2=PointSet(2,RandValue(2));
        %�ӵ㼯��ѡȡ������������ֱ��
        slope=(y2-y1)/(x2-x1);
        intercept=y1-slope*x1;
        %�������е㵽��ֱ�ߵľ���֮��
        for j=1:NonZeroNum+1
            tempSumDistPointsToLine=abs(slope*PointSet(1,j)-PointSet(2,j)+intercept)...
            /sqrt(slope^2+1);
        end
        if tempSumDistPointsToLine<SumDistPointsToLine
            SumDistPointsToLine=tempSumDistPointsToLine;
            FinalSlope=slope;
            FinalIntercept=intercept;
        end
        NumFlag(k)=NonZeroNum;
        Result(1,k)=FinalSlope;
        Result(2,k)=FinalIntercept;
    end
end




function [RandValue]=UserRand(NonZeroNum)
    point1Index=1;
    point2Index=1;
    while point1Index==point2Index||point1Index==0||point2Index==0%��֤���ֱ�ߵ������㲻����ͬ�ĵ㣬��Ȼ�������㱨NaN
        point1Index=round(NonZeroNum*rand(1,1));%roundΪ��������ȡ��
        point2Index=round(NonZeroNum*rand(1,1));
    end
    RandValue=[point1Index point2Index];
end %�Ӻ���

end%������








