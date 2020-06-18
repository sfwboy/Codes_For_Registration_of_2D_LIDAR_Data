clc
clear
[data1]=textread('data1.txt');
[data11]=CoordinateTran1(data1');
[data2]=textread('data2.txt');
[data22]=CoordinateTran1(data2');
figure(1)
hold on
for i=1:size(data11,2)
        plot(data11(1,i),data11(2,i),'g.');   
end
for i=1:size(data22,2)
        plot(data22(1,i),data22(2,i),'r.');   
end
hold off