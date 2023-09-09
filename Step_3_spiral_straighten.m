function picking_file = Spiral_Straighten(imagefile,interpolation_file, mode, round,halfbox, halfWidth, dia, pixA)
% This function is to straighten curved spirals and to calculate much
% accurate picking points. 

% interpolation_file:   is the source interpolation file;
% mode:         1 is for using interpolation data; 2 is for using linspacing data;
% round:        1 is for straightening curved spirals; 2 is for both straightening
%               and calculating new picking points; 
% halfbox:      is the padding size and the boxed image size;2*halfbox must be divisible by 4
% halfWidth:    is the slice size;
% dia:          is the estimated size of spiral diameter to find histogram
%               peak easily;

% halfbox>halfWidth>dia;



roted = [];

ang = [];

im = [];

invert = 1;
% invert=0;

%% find the extension of input file and then open it using right script
m_pth = fileparts(which('Spiral_Straighten.m'));   
[parentdir,~,~] = fileparts(m_pth);   % parent dir
image_path =[parentdir '/Image/'];   %Image path

imagefile_ext = imagefile(strfind(imagefile,'.')+1:end);
    
if strcmp(imagefile_ext,'tif') 
    [Orig_img] = imread([image_path '/' imagefile]); % for linux
    %[Orig_img] = imread([image_path '\' imagefile]);  %for windows
else   % '.mrc' 'jpg'
     [Orig_img s] = ReadEMFile([image_path '/' imagefile]); % for linux
   % [Orig_img s] = ReadEMFile([image_path '\' imagefile]); %for windows
end
if size(Orig_img,3)==3
    Orig_img = rgb2gray(Orig_img);
end;
mn = size(Orig_img,1);
mm = size(Orig_img,2);
% % invert density maps
if (invert == 1) && (round==1)
    a=[mn mm];
    white_mat = single(255.*ones(a));
    Orig_img = white_mat - single(Orig_img);
end
    
%% set paths and file names
pth = ffpath([interpolation_file '.mat']);   
Inter_data =  load([pth '/' interpolation_file '.mat']);

%% choose between Interpolation data or spacing data;
if mode == 1
    slope = Inter_data.slope;
    Points_position = Inter_data.Inter_points;
    seglen = Inter_data.seglen;
else
    slope = Inter_data.linspace_slope;
    Points_position = Inter_data.linspace_points;
    seglen = Inter_data.step_size;
end

diff = zeros(size(Points_position,1)-2,1);

%% pad the original image in case of bigger size of box
im = zeros(mn+2*halfbox,mm+2*halfbox);
im((halfbox+1):(halfbox+mn), (halfbox+1):(halfbox+mm)) = Orig_img;

%% box the image
for i=1 : size(Points_position,1)-1
    cx = Points_position(i,1);
    cy = Points_position(i,2);
    boxim = im(ceil(cy): ceil(cy+2*halfbox-1),ceil(cx): ceil(cx+2*halfbox-1));
    i;

    %% find the rotation angles for each interpolation points
    if Points_position(i,1) >= Points_position((i+1),1)
            ang(i) = 90 -  (-atand(slope(i)));
    else
            ang(i) = 270 - (-atand(slope(i)));   
    end

    %% determine the possible wrong calculation at angles around 0 or 180
    if i>= 2 && abs(ang(i)-ang(i-1))<= 270 && abs(ang(i)-ang(i-1))>= 90
            ang(i) = ang(i) -180;
            if ang(i) <=0
                ang(i)=ang(i)+360;
            end
    end

    ang(i);
    if i >=2
        diff(i-1) = min (abs(ang(i)-ang(i-1)),abs(abs(ang(i)-ang(i-1))-360));
    end

    rot_im = imrotate(boxim,ang(i),'bilinear','crop');

    slice = rot_im((ceil(halfbox-seglen/2)):(ceil(halfbox+seglen/2)-1),halfbox-halfWidth:halfbox+halfWidth);


%     Diff_angle = diff(2:end)';
% 
%     figure(3);
%     plot(Diff_angle(:,1),'k');hold on;

    %% fine shift of each slice
    if round ==1

        for j=1:2*halfWidth
            HistGram(j) = sum(slice(:,j),1);
        end
%         plot(HistGram);

        HistGramNew = HistGram (halfWidth-dia+1:halfWidth+dia);

        Peak = max(HistGramNew);
        PeakY = find (HistGramNew==Peak);
        aa = size (PeakY,1);
        dis = PeakY(ceil(aa/2)) - dia - 1; %dis is the shift distance

        % calculate the new points position in old coordinate system
        Points_position(i,:) = [cx+(dis)*cos(ang(i)*pi/180)  cy+(dis)*sin(ang(i)*pi/180)];    

        % shift slice to new center
        slice = rot_im((ceil(halfbox-seglen/2)):(ceil(halfbox+seglen/2)-1),(halfbox-halfWidth+dis):(halfbox+halfWidth+dis));

    end
    
%% save rotated image
roted = [roted; slice];
    
end

%    Diff_angle = diff(2:end)';
% 
%     figure(3);
%     plot(Diff_angle(:,1),'k');hold on;

interpolation_file = strrep(interpolation_file,'Interpolation','Straighten');

%% save new points positions
%if round ==1
    save ([pth '/' interpolation_file '.mat'],'Points_position','diff','roted');
    picking_file = (interpolation_file);
%end

%% save straighten file into .mrc and .tif files
WriteMRC (roted,pixA,[pth '/' interpolation_file '.mrc'])
imwrite(roted, [pth '/' interpolation_file '.tif'], 'tif')
% save ([pth '/' roted '.mat'],'roted');

%% pack the small area into much larger area to get strong FFT signal
nx = size (roted,1);
ny = size (roted,2);
Pad_im = zeros(nx,ceil(nx/2));
Pad_im(:, ceil((nx/2-ny)/2):ceil((nx/2+ny)/2-1)) = roted;
WriteMRC (Pad_im,pixA,[pth '/' interpolation_file '_pack.mrc']);
imwrite(Pad_im, [pth '/' interpolation_file '_pack.tif'], 'tif');




close all
% figure;SetGrayscale;imagesc(roted);
% imcontrast(gca);

%for histogram
[aa bb]=size(roted);

for i=1:bb
    cc(i) = sum(roted(:,i),1);
end

cc;
% figure(2);
% plot(cc);

cc= cc';

%save ([pth '/diameter.mat'],'cc');
 
