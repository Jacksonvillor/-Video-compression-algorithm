clc;clear all;
fileName = 'cat_video.avi' ;    %traffic.avi  D:666.avi
obj = VideoReader(fileName);
numFrames = obj.NumberOfFrames;% 帧的总数
II=rgb2gray(read(obj,1));
[M,N]=size(II);
File=zeros(M,N,numFrames);
File(:,:,1)=II(:,:);
for oo=1:3
X=input('请输入窗口大小:   ');
Z=input('请输入搜索窗口半径大小:   ');
Y=input('请输入相邻I帧的间隔：');
%WWW=floor(numFrames/Z);
%MOBAN=ones(X,X)/(X^2);  
IZhenWeiZhi=1;
for k = 2:numFrames% 读取数据
     frame =rgb2gray(read(obj,k));
     frameI=rgb2gray(read(obj,IZhenWeiZhi));
     if   ~(mod(k,Y)) %进行桢内预测
         if k+Y<numFrames
         IZhenWeiZhi=k+Y;
         else
            IZhenWeiZhi=k;
         end
         [Izhen,PSNR,CompressionRatio]=IzhenYaSuo(frameI);
         File(:,:,k)=Izhen(:,:); 
         ZPSNR(k)=PSNR;
         YASUOLV(k)=CompressionRatio;    
     else   
         FrameI=frameI;
         [Pzhen,QQQ]=PzhenGuJi(frame,FrameI,X,Z);
         YASUOLV(k)=QQQ;
       %Pzhen=conv2(double(Pzhen),double(MOBAN),'same');%%%%%%%%%%%%%%低通%%%%%%%%%%%%%
       File(:,:,k)=Pzhen(:,:);
       A=rgb2gray(read(obj,k));
       PP=double(A)-double(Pzhen);
       MSE=sqrt(sum(PP(:).^2)/(M*N)); % 均方误差：指参数估计值与参数真值之差平方的期望值，记为MSE
       ZPSNR(k)=10*log(255*255/MSE)/log(10);   
     end;%桢间相似度预测完毕
 end;

 subplot(2,2,oo);
 X=1:numFrames;
 for i=1:numFrames
     Y(i)=ZPSNR(i);
 end;
 save File;
 YaSuoLv=sum(YASUOLV)/(numFrames-1);
plot(X,Y,'r');
xlabel('当前帧数');
ylabel('峰值信噪比');
title({'视频压缩率为：',YaSuoLv});
implay(uint8(File));
%implay('F:MK.avi');
end;