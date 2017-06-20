function [Izhen,PSNR,CompressionRatio]=IzhenYaSuo(frameI)
%%% 开始处理；
%DCT 变换
OriginalImage=double(frameI); % 图像数据类型转换
[am,an]=size(OriginalImage);% 得到图像的大小
ImageSub=OriginalImage-128;
fun1=@dct2;
TCM=blkproc(ImageSub,[8,8],fun1); % 使用dct2函数进行二维DCT变换，得到变换系数矩阵TCM
%3.量化
Q=[16 11 10 16 24 40 51 61
12 12 14 19 26 58 60 55
14 13 16 24 40 57 69 56
14 17 22 29 51 87 80 62
18 22 37 56 68 109 103 77
24 35 55 64 81 101 113 92
49 64 78 87 103 121 120 101
72 92 95 98 112 100 103 99];% 亮度量化表
TCM_Q=blkproc(TCM,[8,8],'round( x./P1)',Q);% 对图像进行量化，得到量化后的TCM，即TCM_Q。
%Z形扫描
TCM_Q_col=im2col(TCM_Q,[8,8],'distinct'); % 将每个8*8 数据块的量化系数排成列向量，得到64* 数据块总数大小的矩阵TCM_Q_col。
Num_col=size(TCM_Q_col,2);% 得到TCM_Q_col的列数，即数据块的个数Num_col。 
order=[1 9 2 3 10 17 25 18 ...
11 4 5 12 19 26 33 41 ...
34 27 20 13 6 7 14 21 ...
28 35 42 49 57 50 43 36 ...
29 22 15 8 16 23 30 37 ...
44 51 58 59 52 45 38 31 ...
24 32 39 46 53 60 61 54 ...
47 40 48 55 62 63 56 64];
TCM_Q_colZ=TCM_Q_col(order,:);% 用z 型扫描方式对变换系数重新排列
%5.编码
%5.1直流编码，dc为直流系数表，dcdpcm为直流差值编码表
dc=zeros(Num_col,1);
dcdpcm=zeros(Num_col,1);
for j=1:Num_col
dc(j)=TCM_Q_colZ(1,j); % 将DC 系数排列到一个矢量中
end
dcdpcm(1)=dc(1);
for j=2:Num_col
dcdpcm(j)=dc(j)-dc(j-1); % 求DC 系数的DPCM 编码
end
dcdmax=max(dcdpcm); %最大直流
dcdmin=min(dcdpcm); %最小直流
dch=histc(dcdpcm,[dcdmin:dcdmax]); %统计各个值的直方图
dcnum=length(dcdpcm);
dcp=dch/dcnum; %计算各个值的概率
dcsymbols=[dcdmin:dcdmax]; %直流分量值
[dcdict,dcavglen]=huffmandict(dcsymbols,dcp); %生成字典dcdict，计算平均码长
dcencoded=huffmanenco(dcdpcm,dcdict); % 对DC 系数的DPCM 进行Huffman 编码，得到直流编码dcencoded

%5.2交流编码
% 将非零AC元素重新排列放到ac中,每一列均以eob 作为结束,共有count个非零元素 
eob=max(ImageSub(:))+1; % 创建一个块结束符号
num=numel(TCM_Q_colZ)+size(TCM_Q_col,2);
ac=zeros(num,1);
count=0;
for j=1:Num_col
i=max(find(TCM_Q_colZ(:,j)));%find 函数为寻找yy 函数中非零元素的位置，max 函数为取里面的最大值，若无非零元素或者为空，返回empty
if isempty(i)
i=1;
end
p=count+1;
q=p+i-1;
if i==1
ac(q)=eob;
end
ac(p:q)=[TCM_Q_colZ(2:i,j);eob];
count=q;
end
ac((count+1):end)=[];% 删除ac中的无用元素
acmax=max(ac); %最大交流
acmin=min(ac); %最小交流
ach=histc(ac,[acmin:acmax]); %统计各个值的直方图
acnum=length(ac);
acp=ach/acnum; %计算各个值的概率
acsymbols=[acmin:acmax]; %交流分量值
[acdict,acavglen]=huffmandict(acsymbols,acp); %生成字典dcdict，计算平均码长
acencoded=huffmanenco(ac,acdict); % 对AC 系数进行Huffman编码，得到交流编码acencoded


%解码
dcdecoded=huffmandeco(dcencoded,dcdict); %直流Huffman解码
%根据直流解码恢复直流分量，并放入TCM_Q_colZ_Rec的第一行
TCM_Q_colZ_Rec(1,1)=dcdecoded(1);
for i=2:Num_col
TCM_Q_colZ_Rec(1,i)=TCM_Q_colZ_Rec(1,i-1)+dcdecoded(i); % 计算第i列直流分量，并将直流分量放入TCM_Q_colZ_Rec的第i列第1行。
end
acdecoded=huffmandeco(acencoded,acdict); %交流Huffman解码
%根据交流解码恢复交流分量，放入TCM_Q_colZ_Rec的第2-64行
j=1; %j用来记录第几列
k=2; %k用来记录第几行
maxk=1;
count=0; %count用来记录连续不等于eob的个数，当count=63时，下一个eob仅作为结束符，不解码。
for i=1:size(acdecoded)
if acdecoded(i)==eob
TCM_Q_colZ_Rec(k:64,j)=0;
j=j+1;
k=2;
else
TCM_Q_colZ_Rec(k,j)=acdecoded(i);
k=k+1;
end
end
%反Z型扫描
order2=[
1 3 4 10 11 21 22 36 ...
2 5 9 12 20 23 35 37 ...
6 8 13 19 24 34 38 49 ...
7 14 18 25 33 39 48 50 ...
15 17 26 32 40 47 51 58 ...
16 27 31 41 46 52 57 59 ...
28 30 42 45 53 56 60 63 ...
29 43 44 54 55 61 62 64];
TCM_Q_col_Rec= TCM_Q_colZ_Rec(order2,:); %用反z 型扫描方式对变换系数重新排列
TCM_Q_Rec=col2im(TCM_Q_col_Rec,[8,8],[am,an],'distinct'); % 将TCM_Q_col_Rec的每个列向量排成8*8数据块，将矩阵变为图像size。
%反量化
TCM_Rec=blkproc(TCM_Q_Rec,[8,8],'round( x.*P1)',Q);
%反DCT变换
fun2=@idct2;
ImageSub_Rec=blkproc(TCM_Rec,[8,8],fun2);% 使用idct2函数进行二维反DCT变换，得到ImageSub_Rec
ReconImage=double(ImageSub_Rec)+128;
Izhen=uint8(ReconImage);

encoded_lenght=numel(dcencoded)+numel(acencoded); %编码长度
AverageBit=encoded_lenght/am/an; %计算编码比特率（每个像素所占的比特数）
CompressionRatio=am*an*8/encoded_lenght; %计算压缩倍数（原图大小与压缩后的比值）
e=double(Izhen)-double(frameI);
MSE=sqrt(sum(e(:).^2)/(an*am)); % 均方误差：指参数估计值与参数真值之差平方的期望值，记为MSE
PSNR=10*log(255*255/MSE)/log(10); %计算峰值信噪比(dB)
end