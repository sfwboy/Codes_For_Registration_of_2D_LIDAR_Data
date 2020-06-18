function [SetOut2]=SplitData1(SetIn,LenThesshold,LineThreshold)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ�1.������ֵThreshold������ĵ㼯SetIn����ֱ�߶ε��жϣ������ֱ��
%           ���ݣ�С����ֵ���򲻴����������ֱ�ߣ����ߣ�������ֵ���������
%           ���ָֻ�ָ�һ�Ρ�
%          2.�ָ����õ�һ���Ӽ����Ӽ������ֱ�߶γ��ȴ���LineThreshold����
%          ΪSetOut�����
%���룺�㼯SetIn��������ֵLenThesshold��ֱ����ֵLineThreshold
%     ���У�1��SetInΪ2��m��n����ά����n����2��m�Ķ�ά���������2Ϊ��ά��
%              ���������mΪ��ά�����������
%              ����ĵ�i��2��m�Ķ�ά�����ʾΪ��
%              SetIn[:,:,i]=|x1 x2 ... xm|
%                           |y1 y2 ... ym|
%              x,y�ֱ�Ϊ�ѿ�������ϵ�еĺ�������
%           2��LenThessholdΪ1��1����ֵ������ֱ�߳�����ֵ
%           3��LineThresholdΪ1��1����ֵ������ֱ�߶��жϵ���ֵ
%��������ݼ�SetOut
%     ���У�1��SetOutΪ2��p��q����ά����q����2��p�Ķ�ά���������2Ϊ��ά��
%              ���������pΪ��ά�����������
%              ����ĵ�j��2��p�Ķ�ά�����ʾΪ��
%              SetOut[:,:,j]=|x1 x2 ... xj|
%                            |y1 y2 ... yj|
%              x,y�ֱ�Ϊ�ѿ�������ϵ�еĺ�������
%���ߣ�Shaofeng Wu 
%ʱ�䣺2018.09.10
%���䣺shaofeng693@126.com
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
    %Step1���жϵ�ǰ�㼯�Ƿ�Ϊֱ�߶μ��ϣ�������ָ�ָ��Ϊ����β�������
    %       ��ֱ����Զ�ĵ�
    Add0=0;
    %*********************************************************************
    %Step1.1����Ͼ�����ǰ���ݼ���β�����ֱ�ߣ�y=ax+b
    %��������ʽ����ֱ�߲���
    lineA=(data1(2,1)-data1(2,len))/(data1(1,1)-data1(1,len));
    lineB=data1(2,1)-data1(1,1)*lineA;
    distMax=0;
    maxFlag=0;
    %*********************************************************************
    %Step1.2�����㵱ǰ�������е㵽ֱ�ߵľ��룬�������ֵ�Ͷ�Ӧ������
    for i=1:len
        %�㵽ֱ�߾��빫ʽ��ע�⣺y=ax+b��Ax+By+C=0֮��Ĳ����任
        distTemp=abs(lineA*data1(1,i)-data1(2,i)+lineB)/sqrt(lineA^2+1);
        if distMax<distTemp
            distMax=distTemp;
            maxFlag=i;
        end   
    end
    tempData=0;
    %*********************************************************************
    %Step1.3���ָ��
    if distMax>LineThreshold%������ֵ�����Ϊ�����Ӽ�
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
    else %С����ֵ���򲻷ָ�
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
    %Step2���жϼ��ϰ���ֱ�߶γ��ȣ�������ֵLenThesshold�����
    if sqrt((SetOut1(1,1,i)-SetOut1(1,add,i))^2+(SetOut1(2,1,i)-SetOut1(2,add,i))^2)>LenThesshold
        for k=1:add
            SetOut2(1,k,Add0)=SetOut1(1,k,i);
            SetOut2(2,k,Add0)=SetOut1(2,k,i);
        end
        Add0=Add0+1;
    end
end





