function [Result_Ls1,NumFlag1]=LeastSquareLine1(SetOut1)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%�������ܣ�����㼯��������С���˷���ϵ㼯�а�����ֱ�߶Σ����ֱ�߲���
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

%��С���˷����ֱ��
[Result_Ls1,NumFlag1]=LeastSquareLine(SetOut1);
%��Ϊ��С���˷����ڼнǽӽ�90�ȵ�ֱ����ϴ��������Խ������²���
%%�ж�����ֱ�ߵ�б�ʣ��������arctan2=63.43���򽫸÷ָ�㼯��ת90�ȣ�
%Ȼ����ʹ����С���˷����ֱ�ߣ�����ٰ�ֱ����ת����
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
%�����ת���ֱ��
[Result_Ls1Tra0,NumFlag1Tra]=LeastSquareLine(Data1);

Line1StartPoint=zeros(2,size(Data1,3));
for j=1:size(Data1,3)
    if abs(Data1(1,1,j)-Data1(1,NumFlag1Tra(j),j))>5
        if TranFlag1(j)~=0
             %ÿ��ֱ��ѡȡ��ʼ�㣬�������������ֱ�߷������
            Line1StartPoint(1,j)=Data1(1,1,j);%��ʼ�������
            Line1StartPoint(2,j)=Result_Ls1Tra0(1,j)*Data1(1,1,j)+Result_Ls1Tra0(2,j);%��ʼ��������
        end
    else
        if TranFlag1(j)~=0
             %ÿ��ֱ��ѡȡ��ʼ�㣬�������������ֱ�߷������
            Line1StartPoint(1,j)=(Data1(2,1,j)-Result_Ls1Tra0(2,j))/Result_Ls1Tra0(1,j);%��ʼ�������
            Line1StartPoint(2,j)=Data1(2,1,j);%��ʼ��������
        end
    end 
end

%�Ѿ�֪����ת��ֱ�ߵ�б��Result_Ls1Tra��Result_Ls2Tra����ת����TranMatrix1
%�Լ�ֱ���ϵ�һ����ʼ��Line1StartPoint��Line2StartPoint
%���Կ��������תǰֱ�ߵ�б��K�ͽؾ�

%��Ϊ��ת90�ȣ�����б��K1=-1/Result_Ls1Tra
for i=1:size(TranFlag1,2)
    if TranFlag1(i)~=0
        Result_Ls1(1,i)=1/Result_Ls1Tra0(1,i);
        temp=Line1StartPoint(:,i)'*(-TranMatrix1);%��ֱ���ϵĵ���ת
        Result_Ls1(2,i)=temp(2)-Result_Ls1(1,i)*temp(1);
    end
end