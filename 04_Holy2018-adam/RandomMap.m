


function [dstImg]=RandomMap(boxNum,boxDst)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%功能：产生由随机排列的方块组成的地图
%输入：方块个数 boxNum
%      方块间的距离 boxDst
%输出：640*480大小的二值图像dstImg
%作者：Shaofeng Wu 
%时间：2018.04.20
%邮箱：shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************


%*********************************************************************
%*********************************************************************
%Step1:创建二值图像，并且在图像上绘制第一个矩形
srcImg=zeros(480,640);   %创建640*480大小的图片srcImg
dstImg=im2bw(srcImg);    %将图片srcImg转为二值图（黑白图）
modelImg=dstImg;         %创建modelImg
left=0;
right=0;
top=0;
down=0;
while (right-left<boxDst||down-top<boxDst)||(right-left>350||down-top>200)%当方块的长宽满足条件是才终止循环
    rectangleX=int32(640*rand());   %随机生成矩形左上角顶点横坐标
    rectangleY=int32(480*rand());   %随机生成矩形左上角顶点纵坐标
    width=int32(400*rand());        %随机生成矩形宽度
    height=int32(200*rand());       %随机生成矩形高度
    left=rectangleX;
    if left<1               %超出图片范围时设为1
        left=1;
    end
    right=rectangleX+width; %超出图片范围时设为640   
    if(right>640)
        right=640;
    end
    top=rectangleY;         %超出图片范围时设为1
    if(top<1)
        top=1;
    end
    down=rectangleY+height;
    if down>480             %超出图片范围时设为480
        down=480;
    end
end  
for i=top:down              %在图片上绘制矩形
    for j=left:right
        dstImg(i,j)=1;
    end
end
%*********************************************************************
%*********************************************************************
%Step2:创建模板图像，用来标记被矩形占据的像素点
[modelImg]=fillModelImg(top,down,left,right,modelImg,boxDst);
%*********************************************************************
%*********************************************************************
%Step3:循环往图像中添加矩形
for i=1:boxNum-1
    [dstImg,top1,down1,left1,right1]=creatBox(dstImg,modelImg);
    %模板
    [modelImg]=fillModelImg(top1,down1,left1,right1,modelImg,boxDst);
end

end%%%%%%%%%%%主函数

%*********************************************************************
%*********************************************************************
%子函数：创建模板图像，用来标记被矩形占据的像素点
%*********************************************************************
%*********************************************************************
function [outImg]=fillModelImg(top,down,left,right,inImg,boxDst)
    outImg=inImg;
    topModel=top-boxDst;
    downModel=down+boxDst;
    leftModel=left-boxDst;
    rightModel=right+boxDst;
    if leftModel<1
        leftModel=1;
    end
    if rightModel>640
        rightModel=640;
    end
    if topModel<1
        topModel=1;
    end
    if downModel>480
        downModel=480;
    end
    for i=topModel:downModel
        for j=leftModel:rightModel
            outImg(i,j)=1;
        end
    end
end%%%%%%%%%%%模板函数

%*********************************************************************
%*********************************************************************
%子函数：创建新的矩形
%*********************************************************************
%*********************************************************************
function [outImg,top1,down1,left1,right1]=creatBox(inImg,modelImg)
    outImg=inImg;
    rectangleX1=int32(640*rand());
    rectangleY1=int32(480*rand());
    while modelImg(rectangleY1,rectangleX1)==1
        rectangleX1=int32(640*rand());
        rectangleY1=int32(480*rand());
    end
    width1=int32(400*rand());
    height1=int32(200*rand());
    if width1<40
        width1=40;
    end
    if height1<20
        height1=20;
    end
    left1=rectangleX1;
    right1=rectangleX1+width1;
    if right1>640
        right1=640;
    end
    while modelImg(rectangleY1,right1)==1
        right1=right1-1;
    end
    top1=rectangleY1;
    down1=rectangleY1+height1;
    if down1>480
        down1=480;
    end
    while modelImg(down1,rectangleX1)==1
        down1=down1-1;
    end
    for i=top1:down1
        for j=left1:right1
            outImg(i,j)=1;
        end
    end

end%%%%%%%[outImg]=creatBox[inImg]

% function [outImg,top,down,left,right]=creatBox(inImg,modelImg)
%     outImg=inImg;
%     left=0;
%     right=0;
%     top=0;
%     down=0;
%     while (right-left<25||down-top<25)||(right-left>350||down-top>200)
%         rectangleX=int32(640*rand());
%         rectangleY=int32(480*rand());
%         width=int32(400*rand());
%         height=int32(200*rand());
%         left=rectangleX;
%         if left<1
%             left=1;
%         end
%         right=rectangleX+width;
%         if(right>640)
%             right=640;
%         end
%         while modelImg(rectangleY,right)==0
%             right=right-1;
%         end
%         top=rectangleY;
%         if(top<1)
%             top=1;
%         end
%         down=rectangleY+height;
%         if down>480
%             down=480;
%         end
%         while modelImg(down,rectangleX)==0
%             down=down-1;
%         end
%     end
%     for i=top:down
%         for j=left:right
%             outImg(i,j)=0;
%         end
%     end
% 
% end%%%%%%%[outImg]=creatBox[inImg]
































