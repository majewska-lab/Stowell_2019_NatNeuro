%% User specified analysis settings menu
% Build user interface for analysis options. Each flag is queried later in
% the script to determine analyses to run and plotting/saving options.
clear all 
close all;
clc

[settings, button] = settingsdlg(...
    'Description'                                       ,'Select Analysis Settings:'                ,...
    'title'                                             ,'Microglial Motility Analysis'             ,...
    'separator'                                         ,'Threshold settings'                       ,...
    {'Enter the timelapse increments (min): ','timeInc'},'5'                                        ,...
    {'Perform threshold test?','threshTestFlag'}        ,[false, false]                              ,...                                                                 ,...
    {'Theshold sensitivity','thresh'}                   ,'0.3'                                      ,...
    'separator'                                         ,'Select soma analyses to run:'             ,...
    {'Count somas','somaCountFlag'}                     ,[false,true]                               ,...
    {'Run soma threshold test','somaThreshFlag'}        ,[false]                                    ,...
    'separator'                                         ,'Select motility analyses to run:'         ,...
    {'Motility index','calcMotFlag'}                    ,[true]                                     ,...
    {'Stability index','calcStabFlag'}                  ,[true]                                     ,...
    {'Instability index','calcInstabFlag'}              ,[true]                                     ,...
    {'Stability histogram','calcHistFlag' }             ,[true]                                     ,...
    'separator'                                         ,'Select pseudopodia analyses to run:'      ,...
    {'Count pseudopodia? ', 'calcPseudoFlag'}           ,[true false]                                ,...
    {'Select minimum pixel area: ','minPix'}            ,'100'                                      ,...
    {'Run pseudopodia timeline?','timePseudoFlag'}      ,[true]                                     ,...
    {'Plot example comparison?','plotExFlag'}           ,[false]                                    ,...
    'separator'                                         ,'Select process coverage analyses to run:' ,...
    {'Calculate process coverage?','areaCovFlag'}       ,[true false]                               ,...
    {'Run process coverage average?','avgCovFlag'}      ,[true]                                     ,...
    {'Run process coverage timeline?','timeCovFlag'}    ,[true]                                     ,...
    'separator'                                         ,'Select output settings:'                  ,...
    {'Save files to other directory','altSaveFlag'}     ,[false]                                    ,...
    {'Save all analysis data?','saveDataFlag'}          ,[true]                                     ,...
    {'Plot all analysis data?','plotDataFlag'}          ,[true false]                               ,...
    {'Save all plots?','savePlotFlag'}                  ,[true]                                     ,...
    'separator'                                         ,[]                                         ...
    );

if strcmp(button,'cancel') || isempty(button)
    fprintf('\n*****************  No settings were entered. Script stopped.  *********************\n\n\n');
    return
end


multicoreFlag = 1; %This flag specifies that multiple cores should be utilized to process the data, if available
close all;
tic
%% Select files to load...

%Opens a user guided file selection window and returns the file name as
%well as the file path. Unless specified otherwise with the "save
%alternative directory," all files saved by this program are saved in the
%path of the file selected here. The file selection window only displays
%tif files.
fprintf('Select tif file...\n');
[imFileName, filePath] = uigetfile('*.tif', 'Select timelapse .tif file to run...');
addpath(filePath);
fprintf('File selected: %s\n',imFileName);
if isempty(imFileName)
    error('No tiff file selected...');
    return;
else
    try
        Img = TIFFStack(imFileName);
        
    catch
        fprintf('Could not find file, select tif file...\n');
        [imFileName, filePath] = uigetfile('*.tif','Select tif file...');
        Img = TIFFStack(imFileName);
    end
end
shortName = strrep(imFileName,'.tif','');

cd(filePath); %Makes the selected files path the current directory
stack = Img(:,:,:); %Creates a working variable copying the orignal tif file
stack=double(stack); %Converts the image stack to double precision.

xpix=size(stack,1); ypix=size(stack,2); zpix=size(stack,3); totpix=xpix*ypix*zpix;
pixArray=reshape(stack,[totpix,1]); subPix=prctile(pixArray,25);
stack=stack-subPix;
stack(stack<0)=0;

settings.rowSize=size(Img,1); %Returns the dimensions of images in the y-axis (rows)
settings.colSize=size(Img,2); %Returns the dimensions of the images in the x-axis (columns)
settings.timePoints=size(Img,3); %Returns the dimension of the images in the z-axis (timepoints)


%% Save alternate directory
%If the user specifies to save the output of this script in another
%directory, this flag triggers and the user picks the new directory
if settings.altSaveFlag==1
    fprintf('You have selected to save files in an alternate directory!\n');
    filePath=uigetdir();
    cd(filePath) %Makes the selected directory the new directory for saving files
    fprintf('Files will be saved in:\n');
    newDir=strcat(pwd,'/','\n');
    fprintf(newDir);
end

%% Identify and count somas
%If the option to count somas was selected this section runs
if settings.somaCountFlag==1
    somaSum=sum(stack,3); %Timepoints are collapsed across the z-axis
    %Normalize the summed image to [0,1]
    Smin = min(somaSum(:));
    Smax = max(somaSum(:));
    nSomaSum = (somaSum - Smin) ./ (Smax - Smin);

    if settings.somaThreshFlag==1
        fprintf('Starting soma threshold test...');
        testThresh=[0.2 0.3 0.4 0.5]; %This array contains default threshold values to test
        testResults=[];
        somaThreshTest=figure('Position',[100 100 1000 1000]);
        hold on
        %The loop runs the thresholding agorithm on the soma projection
        %image with each value in testThresh
        for i=1:length(testThresh)
            [BW,maskedImage] = segmentImage(nSomaSum,testThresh(i),100);
            bwl = bwlabel(BW);
            stats=regionprops(bwl,'area'); % get area of each mode
            A= [stats.Area];
            numSoma=length(A(A>max(A)/4));
            subplot(length(testThresh)/2,length(testThresh)/2,i)
            imagesc(bwl);
            title(['Thresh= ',num2str(testThresh(i)),', Cell # = ',num2str(numSoma)]);
            axis square;
            hold on
        end
        tightfig;
        %This interface allows the user to select the threshold level for
        %soma calculation
        [chooseSoma, button] = settingsdlg(...
        'Description'                                       ,'Select threshold to apply:'    ,...
        'title'                                             ,'Soma threshold test'           ,...
        'separator'                                         ,'Threshold settings'            ,...
        {['Choose: ',num2str(testThresh),' '],'somaThresh'} ,[0]                             ,...
        'separator'                                         ,''                              ...
        );
        settings.somaThresh=chooseSoma.somaThresh;
    else
        testThresh=0.4; %If a threshold test was not selected, the default value is assigned here
        [BW,maskedImage] = segmentImage(nSomaSum,testThresh,100);
        bwl = bwlabel(BW);
        stats=regionprops(bwl,'area');
        A= [stats.Area]; %Produces a list of soma-object areas
        numSoma=length(A(A>max(A)/4)); %Returns only a soma number based on the top 75% areas
        title(['Thresh= ',num2str(testThresh),', Cell # = ',num2str(numSoma)])
        figure;imagesc(bwl)
    end
    data.numSoma=numSoma; %saves the number of somas
    fprintf('...Done!\n');
        
end

%% Binarize stack timepoints

if settings.calcMotFlag==1
    motIndex=zeros(1,settings.timePoints-1); %Creates an empty array based on number of timepoints
    threshStack=[];
    if settings.threshTestFlag==0 %Executes thresholding based on predefined value
      thresh=settings.thresh;
      %Creates a binarized stack based on adaptive, local, pixel-based
      %comparison
      for i=1:settings.timePoints
          binStack=Img(:,:,i);
          binStack=imbinarize(binStack,'adaptive','Sensitivity',thresh);
          threshStack=cat(3,threshStack,binStack);
      end
      if settings.plotDataFlag==1
        figure;
        threshPlot=imagesc(threshStack(:,:,i)); colormap gray; axis square; axis off;
      end
    elseif settings.threshTestFlag==1 %Executes threshold test to allow user to define threshold level
      thresh=[0.2 0.4 0.6 0.8]; motThreshTest=figure('Position',[100 100 3600 450]);
      subplot(1,length(thresh)+1,1);
      imagesc(stack(:,:,1),[-10 100]); axis square; hold on; colormap gray;
      title(['Original image: ']); colorbar;
      for i=1:length(thresh)
          binStack=Img(:,:,1);
          subplot(1,length(thresh)+1,i+1)
          binStack=imbinarize(binStack,'adaptive','Sensitivity',thresh(i));
          imagesc(binStack); colorbar; axis square; colormap gray;
          title(['Thresh= ',num2str(thresh(i)),':'])
      end
      tightfig;
        [chooseThresh, button] = settingsdlg(...
        'Description'                                       ,'Select threshold to apply:'   ,...
        'title'                                             ,'Motility threshold test'      ,...
        'separator'                                         ,'Threshold settings'           ,...
        {['Choose: ',num2str(thresh),' '],'thresh'}         ,[0]                            ,...
        'separator'                                         ,''                              ...
        );
        settings.thresh=chooseThresh.thresh;
      if settings.plotDataFlag==0;
          close all
      end
      for i=1:settings.timePoints
          binStack=Img(:,:,i);
          binStack=imbinarize(binStack,'adaptive','Sensitivity',chooseThresh.thresh);
          threshStack=cat(3,threshStack,binStack);
      end
      figure;
      threshPlot=imagesc(threshStack(:,:,i)); colormap gray; axis square; axis off;
    elseif settings.threshFlag==0
    end
        %Non-threshold motility
        
end
%% Calculate motility index
if settings.calcMotFlag==1;
    %Sets up empty arrays to store pixel counts
    data.extensionPix=[];
    data.retractionPix=[];
    data.stablePix=[];
    data.blankPix=[];
    data.motilityIndex=[];
    %Loops through the stack to count the extension/retraction/stable
    %pixels
    for j=1:size(threshStack,3)-1;
        img1=threshStack(:,:,j); img2=threshStack(:,:,j+1);
        img2(img2>0)=2; %Sets all 1's to 2's in second timepoint
        sumOverlay=img1+img2;%Sums the 1's and 2's to generate an overlay
        red=size(sumOverlay(sumOverlay==1),1);%1's are retraction pixels
        green=size(sumOverlay(sumOverlay==2),1);%2's are extension pixels
        yellow=size(sumOverlay(sumOverlay==3),1);%3's are stable pixels
        blank=size(sumOverlay(sumOverlay==0),1);%0's are blank pixels
        motility=(red+green)/yellow;%Calculate the motility index
        
        %Saves all the results in the data structure
        data.extensionPix=[data.extensionPix green];
        data.retractionPix=[data.retractionPix red];
        data.stablePix=[data.stablePix yellow];
        data.blankPix=[data.blankPix blank];
        data.motilityIndex=[data.motilityIndex motility];
    end
    %Plots the motility index across time if specified by user
    if settings.plotDataFlag==1
        motPlot=figure;
        time=settings.timePoints-1;
        xInc=(1:1:time)*settings.timeInc;
        mot=plot(xInc,data.motilityIndex,'linewidth',3,'color','r');
        %ylim([lowBound, upBound]);
        ylim([0,2]);
        ylabel(['Motility Index']);
        xlabel(['Time (min)']);
        set(gca,'fontname','helvetica','tickdir','out','linewidth',2,'fontsize',16,'fontweight','bold')
        box off;
        title(['Average motility index= ',num2str(mean(data.motilityIndex))])
    end
end

%% Calculate stability index
if settings.calcStabFlag==1;
    %Sets up blank arrays to store data
    data.superstabPix=[];
    stabMat=[];
    %Compares series of 3 images and compares each pixel across time
    for k=1:size(threshStack,3)-2;
        img1=threshStack(:,:,k); img2=threshStack(:,:,k+1); img3=threshStack(:,:,k+2);
        stabTime=zeros(settings.rowSize,settings.colSize);
        for i = 1:settings.rowSize
            for j=1:settings.colSize
                %Saves pixels in a matrix that extended (0->1, and then
                %stayed stable (1->1)
                if img1(i,j)==0 && img2(i,j)==1 && img3(i,j)==1
                    stabTime(i,j)=1;
                end
            end
        end
        %Generates a stack that saves all the super-stable pixels
        stabMat=cat(3,stabMat,stabTime);
    end
    %Checks to see if the motility index has already been calculated
    if settings.calcMotFlag==1;
        extPix=data.extensionPix;
    %Otherwise, calculates the total number of extension pixels here
    elseif settings.calcMotFlag==0;
        extPix=[];
        for j=1:size(threshStack,3)-1;
            img1=threshStack(:,:,j); img2=threshStack(:,:,j+1);
            img2(img2>0)=2;
            sumOverlay=img1+img2;
            green=size(sumOverlay(sumOverlay==2),1);
            extPix=[extPix green];
        end
    end
    data.extensionPix=extPix(1:end-1); %Omits the last extension timepoint because it cannot be compared to a subsequent timepoint
    for h=1:size(stabMat,3)
        tempMat=(stabMat(:,:,h));
        superyellow=size(tempMat(tempMat==1),1); %Calculates the total number of super-stable pixels for each timelapse
        data.superstabPix=[data.superstabPix,superyellow];
    end
    data.stabilityIndex=data.superstabPix./data.extensionPix;%Calculates the stability index by dividing the super-stable pixels by the number of extension pixels
    %Plots the data based on user specified input
    if settings.plotDataFlag==1
        stabPlot=figure;
        time=settings.timePoints-1;
        xInc=(2:1:time)*settings.timeInc;
        stab=plot(xInc,data.stabilityIndex,'linewidth',3,'color','r');
        ylim([0,1]);
        xlim([5,60]);
        ylabel(['Stability Index']);
        xlabel(['Time (min)']);
        set(gca,'fontname','helvetica','tickdir','out','linewidth',2,'fontsize',16,'fontweight','bold')
        box off;
        title(['Average stability index= ',num2str(mean(data.stabilityIndex))])
    end
end

%% Calculate instability index
if settings.calcInstabFlag==1;
    %Sets up blank arrays to store data
    data.instabilityPix=[];
    data.instabIndex=[];
    instMat=[];
    %Compares a series of three images across time
    for k=1:size(threshStack,3)-2;
        img1=threshStack(:,:,k); img2=threshStack(:,:,k+1); img3=threshStack(:,:,k+2);
        instTime=zeros(settings.rowSize,settings.colSize);
        for i = 1:settings.rowSize
            for j=1:settings.colSize
                if img1(i,j)==1 && img2(i,j)==1 && img3(i,j)==0
                    instTime(i,j)=1;
                end
            end
        end
        instMat=cat(3,instMat,instTime);
    end
    if settings.calcMotFlag==1;
        stabPix=data.stablePix;
    elseif settings.calcMotFlag==0;
        stabPix=[];
        for j=1:size(threshStack,3)-1;
            img1=threshStack(:,:,j); img2=threshStack(:,:,j+1);
            img2(img2>0)=2;
            sumOverlay=img1+img2;
            yellow=size(sumOverlay(sumOverlay==3),1);;
            stabPix=[stabPix yellow];
        end
    end
    data.stabilityPix=stabPix(1:end-1);
    for h=1:size(instMat,3)
        tempMat=(instMat(:,:,h));
        instPix=size(tempMat(tempMat==1),1);
        data.instabilityPix=[data.instabilityPix,instPix];
    end
    data.instabIndex=data.instabilityPix./data.stabilityPix;
    if settings.plotDataFlag==1
        instPlot=figure;
        time=settings.timePoints-1;
        xInc=(2:1:time)*settings.timeInc;
        inst=plot(xInc,data.instabIndex,'linewidth',3,'color','r');
        %ylim([lowBound, upBound]);
        ylim([0,0.5]);
        xlim([5,60]);
        ylabel(['Instability Index']);
        xlabel(['Time (min)']);
        set(gca,'fontname','helvetica','tickdir','out','linewidth',2,'fontsize',16, 'fontweight','bold')
        box off;
        title(['Average instability index= ',num2str(mean(data.instabIndex))])
    end
end

%% Calculate stability histogram
if settings.calcHistFlag==1;
    histStack=threshStack;
    for k=2:size(histStack,3)-1;
        for i=1:settings.rowSize
            for j=1:settings.colSize
                if histStack(i,j,k)>=1 && histStack(i,j,k+1)>=1
                    histStack(i,j,k+1)=histStack(i,j,k-1)+histStack(i,j,k);
                end
            end
        end
    end
    maxPix=fibonacci(settings.timePoints);
    for i=1:settings.rowSize
        for j=1:settings.colSize
            if histStack(i,j,settings.timePoints)==maxPix;
                histStack(i,j,:)=0;
            end
        end
    end
    sumStack=sum(histStack,3);
    histVec=reshape(sumStack,1,settings.rowSize*settings.colSize);
    histVec(histVec==0)=[]; histVec(histVec>30)=[];
    allPix=size(histVec,2);
    data.lowStab=size(histVec(histVec<=4),2)/allPix;histVec(histVec<=4)=[];
    data.midStab=size(histVec(histVec<=20),2)/allPix;histVec(histVec<=20)=[];
    data.hiStab=size(histVec(histVec<=88),2)/allPix;histVec(histVec<=88)=[];
    countPix=[data.lowStab data.midStab data.hiStab];
    if settings.plotDataFlag==1
        xtitles={'Low';'Medium';'High'};
        figure; histPlot=bar(countPix);
        axis square;
        set(gca,'fontname','helvetica','tickdir','out','linewidth',3,'fontsize',16,'xticklabel',xtitles,'fontweight','bold')
        set(histPlot,'linewidth',2,'FaceColor','r');
        title('Stability histogram');
        ylabel(['Pixel proportion']); ylim([0,1]);
        xlabel(['Pixel stability']);
    end
end

%% Pseudopodia analysis
if settings.calcPseudoFlag==1
    minPix=settings.minPix;
    extStack=[];
    retStack=[];
    extCount=[];
    retCount=[];
    for j=1:size(threshStack,3)-1;
        img1=threshStack(:,:,j); img2=threshStack(:,:,j+1);
        img2(img2>0)=2;
        diffOverlay=img2-img1;
        extOverlay=diffOverlay; extOverlay(extOverlay<2)=0; extOverlay(extOverlay==2)=1;
        retOverlay=diffOverlay; retOverlay(retOverlay>-1)=0; retOverlay(retOverlay==-1)=1;
        extStack=cat(3,extStack, extOverlay);
        retStack=cat(3,retStack, retOverlay);
    end
    for i=1:size(threshStack,3)-1;
        testThresh=0.5;
        [ext,maskedImageExt] = segmentImage(extStack(:,:,i),testThresh,8);
        [ret,maskedImageRet] = segmentImage(retStack(:,:,i),testThresh,8);
        extl = bwlabel(ext);
        retl = bwlabel(ret);
        statExt=regionprops(extl,'area');
        statRet=regionprops(retl,'area');
        extA=[statExt.Area];
        retA=[statRet.Area];
        extCount=cat(2,extCount,(size(extA(extA>=minPix),2)));
        retCount=cat(2,retCount,(size(retA(retA>=minPix),2)));
    end
    if settings.plotDataFlag==1 && settings.timePseudoFlag==1;
        pseudoPlot=figure;
        time=settings.timePoints-1;
        xInc=(1:1:time)*settings.timeInc;
        extP=plot(xInc,extCount,'linewidth',3,'color','g');
        hold on
        retP=plot(xInc,retCount,'linewidth',3,'color','r');
        ylim([0 100]);
        ylabel(['Pseudopodia counts']);
        xlabel(['Time (min)']);
        set(gca,'fontname','helvetica','tickdir','out','linewidth',2,'fontsize',16,'fontweight','bold')
        box off;
        title(['Extension/Retraction pseudopodia timeline'])
        legend('extension','retraction');
    end
    
    data.pExtTime=extCount;
    data.pRetTime=retCount;
    data.pExtRatio=extCount./(extCount+retCount);
    data.pRetRatio=retCount./(extCount+retCount);
    
     if settings.plotExFlag==1
        examplePlot=figure('Position',[10 10 1200 600]); hold on
        subplot(2,3,[1,4])
        imagesc(diffOverlay); colorbar; axis square; axis off;
        subplot(2,3,2)
        imagesc(extOverlay); axis square; axis off;
        subplot(2,3,3)
        imagesc(retOverlay); axis square; axis off;
        subplot(2,3,5)
        imagesc(extl); axis square; axis off;
        subplot(2,3,6)
        imagesc(retl); axis square; axis off;
        tightfig;
     end
     
end

%% Process area coverage analysis
if settings.areaCovFlag==1
    timeCov=[];
    areaStack=sum(threshStack,3);
    totalPix=settings.rowSize*settings.colSize;
    for i=1:size(threshStack,3)
        pixCov=threshStack(:,:,i);
        pixCov=size(pixCov(pixCov==1),1);
        timeCov=[timeCov pixCov];
    end
    timeCov=timeCov/totalPix;
    coverage=size(areaStack(areaStack>0),1);
    pixCovRatio=coverage/totalPix;
    
    if settings.avgCovFlag==1
        data.totCov=pixCovRatio;
        data.avgCov=mean(timeCov);
    end
    data.timeCov=timeCov
    if settings.timeCovFlag==1
        if settings.plotDataFlag==1
            Cov=figure('Position',[100 100 1600 400]); 
            subplot(1,2,1);
            imagesc(areaStack); axis square; colormap hot; colorbar;
            set(gca,'fontname','helvetica','fontsize',16,'fontweight','bold')
            title(['Process area coverage: ' mat2str(pixCovRatio)]);
            axis off;

            covSub=subplot(1,2,2);
            time=settings.timePoints;
            xInc=(1:1:time)*settings.timeInc;
            cover=plot(xInc,data.timeCov,'linewidth',3,'color','r');
            ylim([0,0.5]);
            ylabel(['Process coverage ratio']);
            xlabel(['Time (min)']);
            set(gca,'fontname','helvetica','tickdir','out','linewidth',2,'fontsize',16,'fontweight','bold')
            box off;
        end
    end
end

%% Saving files and figures
fn=imFileName(1,1:end-4);

if settings.saveDataFlag==1
    save(strcat(fn,'-data.mat'),'data');
    fileID=fopen(strcat(fn,'-data.txt'),'wt');
    outputData=struct2cell(data);
    outputHead=fieldnames(data);
    output=cat(2,outputHead,outputData);
    for h=1:length(outputData)
        outputData{h,1}=outputData{h,1}';
    end
    outputData=outputData';
    for k=1:size(output,1)
        format='%-s\t';
        for j=1:length(output{k,2})
            format=strcat(format,' %-1.2f');
        end
        format=strcat(format,' \n');
        fprintf(fileID,format,output{k,:});
    end;
    fclose(fileID);
    parameters=struct2cell(settings);
    paramHead=fieldnames(settings);
    parameters=cat(2,paramHead,parameters);
    
    save(strcat(fn,'-settings.mat'),'settings');
    fileID2=fopen(strcat(fn,'-settings.txt'),'wt');
    for k=1:size(parameters,1)
        format='%s\t %1.1f \n';
        fprintf(fileID2,format,parameters{k,:});
    end;
    fclose(fileID2);
end

if settings.plotDataFlag==1 && settings.savePlotFlag==1 
    if settings.somaCountFlag==1 & settings.somaThreshFlag==1
        saveas(somaThreshTest,strcat(fn,'-somaThreshTest.tif'));
    end
    if settings.threshTestFlag==1;
        saveas(motThreshTest,strcat(fn,'-thTest.tif'));
    end
    if settings.calcMotFlag==1;
        saveas(motPlot,strcat(fn,'-mIndex.tif'));
    end
    if settings.calcStabFlag==1;
        saveas(stabPlot,strcat(fn,'-sIndex.tif'));
    end
    if settings.calcInstabFlag==1;
        saveas(instPlot,strcat(fn,'-iIndex.tif'));
    end
    if settings.calcHistFlag==1;
        saveas(histPlot,strcat(fn,'-sHist.tif'));
    end
    if settings.calcPseudoFlag==1 && settings.timePseudoFlag==1; 
        saveas(pseudoPlot,strcat(fn,'-pTime.tif'));
    end
    if settings.plotExFlag==1;
        saveas(examplePlot,strcat(fn,'-pEx.tif'));
    end
    if settings.areaCovFlag==1;
        saveas(Cov,strcat(fn,'-covTime.tif'));
    end
end

%close all
%clear all



