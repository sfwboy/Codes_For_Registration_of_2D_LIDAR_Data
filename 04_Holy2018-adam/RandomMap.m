


function [dstImg]=RandomMap(boxNum,boxDst)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%���ܣ�������������еķ�����ɵĵ�ͼ
%���룺������� boxNum
%      �����ľ��� boxDst
%�����640*480��С�Ķ�ֵͼ��dstImg
%���ߣ�Shaofeng Wu 
%ʱ�䣺2018.04.20
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************


%*********************************************************************
%*********************************************************************
%Step1:������ֵͼ�񣬲�����ͼ���ϻ��Ƶ�һ������
srcImg=zeros(480,640);   %����640*480��С��ͼƬsrcImg
dstImg=im2bw(srcImg);    %��ͼƬsrcImgתΪ��ֵͼ���ڰ�ͼ��
modelImg=dstImg;         %����modelImg
left=0;
right=0;
top=0;
down=0;
while (right-left<boxDst||down-top<boxDst)||(right-left>350||down-top>200)%������ĳ������������ǲ���ֹѭ��
    rectangleX=int32(640*rand());   %������ɾ������ϽǶ��������
    rectangleY=int32(480*rand());   %������ɾ������ϽǶ���������
    width=int32(400*rand());        %������ɾ��ο��
    height=int32(200*rand());       %������ɾ��θ߶�
    left=rectangleX;
    if left<1               %����ͼƬ��Χʱ��Ϊ1
        left=1;
    end
    right=rectangleX+width; %����ͼƬ��Χʱ��Ϊ640   
    if(right>640)
        right=640;
    end
    top=rectangleY;         %����ͼƬ��Χʱ��Ϊ1
    if(top<1)
        top=1;
    end
    down=rectangleY+height;
    if down>480             %����ͼƬ��Χʱ��Ϊ480
        down=480;
    end
end  
for i=top:down              %��ͼƬ�ϻ��ƾ���
    for j=left:right
        dstImg(i,j)=1;
    end
end
%*********************************************************************
%*********************************************************************
%Step2:����ģ��ͼ��������Ǳ�����ռ�ݵ����ص�
[modelImg]=fillModelImg(top,down,left,right,modelImg,boxDst);
%*********************************************************************
%*********************************************************************
%Step3:ѭ����ͼ������Ӿ���
for i=1:boxNum-1
    [dstImg,top1,down1,left1,right1]=creatBox(dstImg,modelImg);
    %ģ��
    [modelImg]=fillModelImg(top1,down1,left1,right1,modelImg,boxDst);
end

end%%%%%%%%%%%������

%*********************************************************************
%*********************************************************************
%�Ӻ���������ģ��ͼ��������Ǳ�����ռ�ݵ����ص�
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
end%%%%%%%%%%%ģ�庯��

%*********************************************************************
%*********************************************************************
%�Ӻ����������µľ���
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
































