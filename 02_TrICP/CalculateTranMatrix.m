

function [R,T,theta]=CalculateTranMatrix(set)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：计算两个点集P,X之间的旋转矩阵（R）和平移矩阵（T）
%输入：包含两个一一对应的点集对的矩阵set
%       1.set为4×n矩阵，set(:,i)=| xP  |    
%                                 | yP  |     
%                                 | xX  |    
%                                 | yX  |     
%输出:两个矩阵 旋转矩阵（R）和平移矩阵（T）
%其中： 1.R为2×2矩阵，R=|cos(alpha)   sin(alpha)|
%                        |-sin(alpha)  cos(alpha)|
%       2.T为2×1矩阵，T=|△x|
%                        |△y|
%其中(xP,yP)为待匹配点集P，(xX,yX)模板点集X
%set点集数据基于笛卡尔坐标系
%变量单位均为米制（m）
%作者：Shaofeng Wu 
%时间：2017.11.22
%邮箱：shaofeng693@126.com

%*********************************************************************
%*********************************************************************
%Step1:计算两组点集的中心
pSetCenterPoint=[0,0];%待匹配点集中心（xPCenter,yPCenter）
xSetCenterPoint=[0,0];%模板点集中心(xXCenterX,yXCenter)
xP=set(1,:);          %从4×1维矩阵中读取行向量，便于计算中心
yP=set(2,:);
xX=set(3,:);
yX=set(4,:);
pSetCenterPoint(1)=sum(xP)/numel(xP);   %待匹配点集P中心横坐标
pSetCenterPoint(2)=sum(yP)/numel(yP);   %待匹配点集P中心纵坐标
xSetCenterPoint(1)=sum(xX)/numel(xX);   %模板点集X中心横坐标
xSetCenterPoint(2)=sum(yX)/numel(yX);   %模板点集X中心纵坐标
%*********************************************************************
%*********************************************************************
%Step2:分别计算两个点集P、X每个点到其点集中心的偏移
xPOffset=xP- pSetCenterPoint(1);
yPOffset=yP- pSetCenterPoint(2);
xXOffset=xX- xSetCenterPoint(1);
yXOffset=yX- xSetCenterPoint(2);
%*********************************************************************
%*********************************************************************
%Step3:计算变换矩阵R,T
Numerator=sum(xXOffset.*yPOffset-yXOffset.*xPOffset);   %分子
Denominator=sum(xXOffset.*xPOffset+yXOffset.*yPOffset); %分母
alpha=atan(Numerator/Denominator);                      %偏移角度角度
R=[cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];       %旋转矩阵
T=(xSetCenterPoint')-R*(pSetCenterPoint');              %平移矩阵 
theta=alpha;
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 