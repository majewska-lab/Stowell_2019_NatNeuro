% Vector analysis:
% Created by Rianne D Stowell and Brendan S Whitelaw, Majewska Lab


%Select input folder 
fprintf('Select input folder');
input_dir=uigetdir();
addpath(input_dir);
cd(input_dir);
%Get file list, and remove trailing spaces
input_names=list_directory(input_dir,'.avi');
input_names=strtrim(input_names);
%Run 'goatland' ablation analysis for each file. Results are stored in
%'results' cell array. 'goatland' is the name of the function that runs the
%vector analysis for each movie file. The output of this function is a
%cell array, with each output in a single row. In the 'results' array, each
%row then consists of data from a single video file. 

results=cell(size(input_names,1),8);

for i=1:size(input_names,1)
    vector_output=Vector_fct(input_names{i});
    for j=1:8
        results{i,j}=vector_output{1,j};
    end
    fprintf(strcat(num2str(i),'of',num2str(size(input_names,1))))
end

fprintf('All DONE');

function [goat_milk] = Vector_fct(imFileName)

if isempty (imFileName)
    error ('No file selected...');
end
%The following defines my opened file as the variable Goat%
    Goat = imFileName;
    %The following code creates an optic flow video using the Farneback
    %method which tracks movement over my whole video file%
   vidReader = VideoReader(Goat);
   %set up optical flow object to do the estimate%
   opticFlow = opticalFlowFarneback; 
   %Read in video frames and estimate optical flow of each frame. Display
   %the video frames with flow vectors%
   
   xvel_export=[];yvel_export=[];mag_export=[];orient_export=[];
   %Initialize an empty matrix to load values into.
     
   while hasFrame(vidReader);
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    %Generate array for estimated optic flow
    flow = estimateFlow(opticFlow,frameGray);
    %Pull out velocity values from flow array
    x_vel = [flow.Vx];
    y_vel = [flow.Vy];
    
    xvel_export(end+1,:,:)=x_vel;
    yvel_export(end+1,:,:)=y_vel;
    %Creation of 2-3D matrices consisting of each iteration of data 
    %calculated during run (or for all frames)
    
    magnitudes = flow.Magnitude;
    orientations = flow.Orientation;
    
    mag_export(end+1,:,:)=[;magnitudes];
    orient_export(end+1,:,:)=[;orientations];

   end
%% Relative velocity towards center
%Compute dot product between velocity vector and vector towards center for
%points within a donut-shaped ROI around ablation core
close all
%Define image parameters
n_rows=size(xvel_export,2);
n_cols=size(xvel_export,3);
n_time=size(xvel_export,1);

%Generate ROI mask in shape of donut. 
%1. Draw ellipse and position around ablation core.
video=VideoReader(Goat);
vid1=readFrame(video);
imshow(vid1)
core=impoly; %generates 'polygon' object
pos1=wait(core);
core_mask=createMask(core); %binary mask of core
core_pos=getPosition(core); %gives position of polygon vertices 
%as n by 2 matrix with x,y coordinates... so n_col by n_row
%Calculate core_center in row by column coordinates
core_center_col=(max(core_pos(:,1))+min(core_pos(:,1)))/2;
core_center_row=(max(core_pos(:,2))+min(core_pos(:,2)))/2;
%Calculate core size, similar to 'diameter' of core
core_xsize=max(core_pos(:,1))-min(core_pos(:,1));
core_ysize=max(core_pos(:,2))-min(core_pos(:,2));
core_size=(core_xsize+core_ysize)/2;

% Generate mask that excludes core
mask=ones(n_rows,n_cols);
mask=mask-core_mask;
%Remove core from analysis of velocity vectors
xvel_masked=xvel_export;
yvel_masked=yvel_export;
for i=1:n_time
    xvel_masked(i,:,:)=squeeze(xvel_export(i,:,:)).*mask;
    yvel_masked(i,:,:)=squeeze(yvel_export(i,:,:)).*mask;
end

% Generated normalized relative position matrix for each pixel of the
% Make sure to do this in n_rows by n_cols notation (i.e. y by x)
%Generate matrix with row (y) values relative to core center
row_lin=linspace(1,n_rows,n_rows);
row_lin=transpose(row_lin);
row_mat=repmat(row_lin,1,n_cols);
row_mat_rel=-(row_mat-core_center_row);
%Generate matrix with column (x) values relative to core center
col_lin=linspace(1,n_cols,n_cols);
col_mat=repmat(col_lin,n_rows,1);
col_mat_rel=-(col_mat-core_center_col);
%Concatenate positions matrices (x then y) and normalize relative position
%vectors
rel_pos=cat(3,col_mat_rel,row_mat_rel);
norm_rel_pos=vecnorm(rel_pos,2,3);
norm_row_mat=row_mat_rel./norm_rel_pos;
norm_col_mat=col_mat_rel./norm_rel_pos;
norm_rel_pos_mat=cat(3,norm_col_mat,norm_row_mat);
norm_rel_pos_mat(isnan(norm_rel_pos_mat))=0; %set NaN at origin to zero
%Concatenate velocity matrices
velocity_mat=cat(4,xvel_masked,yvel_masked);
%For each time point, for each pixel, calculate dot product of velocity
%vector and normalized relative position vector
rel_vel_mat=zeros(n_rows,n_cols,n_time); 
for i=1:n_time
    vel=squeeze(velocity_mat(i,:,:,:));
    rel_vel=dot(vel,norm_rel_pos_mat,3);
    rel_vel_mat(:,:,i)=rel_vel;
end

%RIANNE'S Full image analysis, produces outputs  in goat_milk%



%filter steps for positive values over thresh%
sort=(rel_vel_mat)<5;
Bsort=sort-1;
dotsort=Bsort.*(rel_vel_mat);
Absolutegoats=abs(dotsort);

%make histogram%

    edges=linspace(5,40,35);
    N=hist(Absolutegoats,2,edges);

%Getting mean of vectors at each tp%
rowmean=sum(Absolutegoats,2)./sum(Absolutegoats~=0,2);
squeeze(rowmean);
Goatking=nanmean(rowmean);
squeeze(Goatking);

%delete first tp%
Goatking(1)=[];
Goatking(1)=[];
%AOC of mean vectors at each tp%
Goatrider=trapz(Goatking);
%Sum of all the vectors%
Goatarmy=sum(Absolutegoats);
Goatnation=sum(Goatarmy);
img_area=n_rows*n_cols;
Goatnations=Goatnation/img_area;
squeeze(Goatnations);
Goatnations(1)=[];
Goatnations(1)=[];
squeeze(Goatnations);
Goatlord=trapz(Goatnations);

%export values into a cell array
goat_milk=cell(1,8);
goat_milk{1,1}=xvel_export;%array containing x velocities
goat_milk{1,2}=yvel_export;%array containing y velocities
goat_milk{1,3}=rel_vel_mat;%array the relative velocity towards the center of the ROI
goat_milk{1,4}=Goatking;%avg at each tp for all relvel sorted vectors%
goat_milk{1,5}=Goatnations; %sum of all vectors at each tp%
goat_milk{1,6}=Goatrider;%AOC of values ommitting first tp%
goat_milk{1,7}=Goatlord;%AOC of sum values%
goat_milk{1,8}=N; %Histogram of values%
end




