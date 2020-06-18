function [SetOut]=SplitData2(SetIn,LenThesshold,LineThreshold)
%**************************************************************************
%**************************************************************************
%**************************************************************************
%�������ܣ�1.������ֵThreshold������ĵ㼯SetIn����ֱ�߶ε��жϣ������ֱ��
%           ���ݣ�С����ֵ���򲻴����������ֱ�ߣ����ߣ�������ֵ���������
%           ���ָֱ���ָ���Ӽ�С����ֵ����Ϊ��ֱ�ߣ�Ϊֹ��
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
%ʱ�䣺2018.10.19
%���䣺shaofeng693@126.com
%**************************************************************************
%**************************************************************************
%**************************************************************************

Flag=true;
NumSetOut=1;
for i=1:size(SetIn,3)
    data=0;%����
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
    %Step1���жϵ�ǰ�㼯�Ƿ�Ϊֱ�߶μ��ϣ�������ָ�ָ��Ϊ����β�������
    %       ��ֱ����Զ�ĵ�
    while Flag
        for add1=1:size(data,3)
            %*********************************************************************
            %Step1.1����Ͼ�����ǰ���ݼ���β�����ֱ�ߣ�y=ax+b
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
            %��������ʽ����ֱ�߲���
            lineA=(temp(2,1,ii)-temp(2,len,ii))/(temp(1,1,ii)-temp(1,len,ii));%б��
            lineB=temp(2,1,ii)-temp(1,1,ii)*lineA;                            %�ؾ�
            distMax=0;
            maxFlag=0;
            %*********************************************************************
            %Step1.2�����㵱ǰ�������е㵽ֱ�ߵľ��룬�������ֵ�Ͷ�Ӧ������
            for k=1:len
                %�㵽ֱ�߾��빫ʽ��ע�⣺y=ax+b��Ax+By+C=0֮��Ĳ����任
                distTemp=abs(lineA*temp(1,k,ii)-temp(2,k,ii)+lineB)/sqrt(lineA^2+1);
                if distMax<distTemp
                    distMax=distTemp;
                    maxFlag=k;
                end   
            end
            %*********************************************************************
            %Step1.3���ָ��
            if distMax>LineThreshold%������ֵ�����Ϊ�����Ӽ�
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
            else      %С����ֵ�򲻷ָ�
                for jj=1:len
                    data(1,jj,ii)=temp(1,jj,ii);
                    data(2,jj,ii)=temp(2,jj,ii); 
                end
                Add=ii;
            end
            end
        end
        %ѭ����ֹ��������ǰ���ϲ����ٷָ������Ϊtemp��data�����а����Ķ�ά���������һ��
    if size(temp,3)~=size(data,3)
        Flag=true;
    else
        Flag=false;
    end
    end
    %*********************************************************************
    %*********************************************************************
    %Step2���жϼ��ϰ���ֱ�߶γ��ȣ�������ֵLenThesshold�����
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




