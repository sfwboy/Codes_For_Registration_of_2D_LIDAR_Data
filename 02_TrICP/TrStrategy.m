function [set]=TrStrategy(DataIn,RateNp)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%���ܣ� TrICP�㷨���Ե�ʵ��
%���ߣ�Shaofeng Wu 
%ʱ�䣺2019.04.24
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************
% DataIn=[2 4 6 1 -2 0 3 9 49 0 8 3 4 2 3;4 5 2 1 0 31 45 8 0 17 33 44 32 24 56; ...
%          23 45 11 42 7 8 98 23 45 23 45 46 8 4 2;-1 24 4 2 212 45 6 89 98 3 2 5 7 -14 -10 ];
% RateNp=0.7;
 
DistSet=zeros(2,size(DataIn,2));%DistSet(1,i)Ϊ��i����Ӧ�����Եľ��룬DistSet(2,i)Ϊ�±�����
%����ÿ����Ӧ��Ե�ŷʽ����
for i=1:size(DataIn,2)
    DistSet(1,i)=(DataIn(1,i)-DataIn(3,i))^2+(DataIn(2,i)-DataIn(4,i))^2;
    DistSet(2,i)=i;
end%for i=1:size(DataIn,2)
%ð������,����
% DistSet=[2 4 6 1 -2 0 3 9;1 2 3 4 5 6 7 8];
for i=1:size(DistSet,2)
    for j=1:size(DistSet,2)-i
        if DistSet(1,j)>=DistSet(1,j+1)
            temp1=DistSet(1,j);
            DistSet(1,j)=DistSet(1,j+1);
            DistSet(1,j+1)=temp1;
            
            temp2=DistSet(2,j);
            DistSet(2,j)=DistSet(2,j+1);
            DistSet(2,j+1)=temp2;
        end
    end%for j=i:size(DistSet,2)-i-1
end%i=1:size(DistSet,2)

%ǰN��������,N=Total*RateNp(�ܵĵ���������ı���)
 N=floor(size(DataIn,2)*RateNp);%����ȡ��
 set=zeros(4,N);
 for i=1:N
     set(1,i)=DataIn(1,DistSet(2,i));
     set(2,i)=DataIn(2,DistSet(2,i));
     set(3,i)=DataIn(3,DistSet(2,i));
     set(4,i)=DataIn(4,DistSet(2,i));
 end




