clear all;
left_im=imread('C:\Users\Melike Nur Mermer\Desktop\Odev2\viewL.png');
right_im=imread('C:\Users\Melike Nur Mermer\Desktop\Odev2\viewR.png');
left_im=rgb2gray(left_im);
right_im=rgb2gray(right_im);
left_d=im2double(left_im);
right_d=im2double(right_im);
% [yuk,gen]=size(left_im);
% disparityRange = [-6 10];
% disparityMap = disparity(left_d,right_d,'BlockSize',...
%     5,'DisparityRange',disparityRange);
% disparityMap1=disparityMap/max(max(disparityMap));
% imshow(disparityMap1);
% c = disparity(left_im,right_im);
w=8;
d1=75;
% left_d=padarray(left_d, [0 d1]);%zero-padding
% right_d=padarray(right_d, [0 d1]);
[yuk,gen]=size(left_d);

% disp=zeros(size(left_d));
% maks=zeros(size(left_d));
tic
for i=1:yuk
    minr = max(1, i-w);
    maxr = min(yuk, i+w);
for j=1:gen
% p=left_d(i,j);
minc = max(1, j-w);
maxc = min(gen, j+w);

		mind = 0; %düþey yönde kayma yok
        maxd = min(d1, gen-maxc);
        
maxcorr=0;
normcorrel=0;
norm_payda=0;
% corr=0;
corr_mat=right_d(minr:maxr, minc:maxc);
numBlocks = maxd - mind + 1;
blockDiffs = zeros(numBlocks, 1);
correl=zeros(numBlocks,1);
ort_corr=mean(mean(corr_mat));
norm1=sum(sum((corr_mat-ort_corr).^2));
for d=mind:maxd
wind_mat=left_d(minr:maxr, (d+minc):(d+maxc));
blockIndex = d - mind + 1;

% blockDiffs(blockIndex, 1) = sum(sum((corr_mat-wind_mat).^2));
blockDiffs(blockIndex, 1) = corr2(wind_mat,corr_mat);
% corr=sum(sum(corr_mat-wind_mat).^2);

ort_wind=mean(mean(wind_mat));
norm2=sum(sum((wind_mat-ort_wind).^2));
norm3=sqrt(norm1*norm2);
normed(blockIndex)=blockDiffs(blockIndex, 1)./norm1;
%     for k=1:2*w+1
%     for l=1:2*w+1
% %         correl(blockIndex)=correl(blockIndex)+left_d(i+k,j+l)*right_d(i+k,j+l-d));
%         normcorrel(blockIndex)=normcorrel(blockIndex)+(left_d(i+k,j+l)-ort_wind)*(right_d(i+k,j+l+d)-ort_corr); 
%         norm_payda(blockIndex)=norm_payda(blockIndex)+((left_d(i+k,j+l)-ort_wind).^2)*((right_d(i+k,j+l+d)-ort_corr).^2);        
%     end
%     end
% normed(blockIndex)=normcorrel(blockIndex)/sqrt(norm_payda(blockIndex));
end
% [temp, sortedIndeces] = sort(blockDiffs,'descend');
[temp, sortedIndeces] = sort(normed,'descend');
bestMatchIndex = sortedIndeces(1, 1);
bm = bestMatchIndex + mind -1;
DbasicSubpixel(i, j) = bm;
% [maks(i,j) disp(i,j)]=max(CC);
% CC(i,j)=disp;
end
end
DbasicSubpixel=uint8(DbasicSubpixel);
imshow(DbasicSubpixel);
toc

%error
load('disp.mat');
figure;
imshow(L);
A_error = double(L)-double(DbasicSubpixel);
min_err=min(min(A_error));
max_err=max(max(A_error));
ort_err=mean(mean(A_error));
std_err=std2(A_error);

% disp=uint8(disp*255);
% CC=uint8(CC);
% imshow(CC);
% fark=left_d-right_d;
% ssd=sum(fark(:).^2);%sum of squared distances
% 
% c=left_d.*right_d;
% cc=sum(sum(c));%cross correlation
% 
% n=left_d.*left_d;
% norm=sum(sum(sum(n)));%for normalization
% 
% nc=cc/norm;%normalized cross correlation