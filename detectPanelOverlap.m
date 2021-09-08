%% Detect Panel Overlap, following twitter feeds:
%
originalFeeds={'https://twitter.com/MicrobiomDigest/status/1433512130406408194';...
    'https://twitter.com/MicrobiomDigest/status/1416115648648908817';...
    'https://twitter.com/MicrobiomDigest/status/1416771640411791363';...
    'https://twitter.com/MicrobiomDigest/status/1416774642069610500';...
    'https://twitter.com/MicrobiomDigest/status/1417391696149635084';...
    'https://twitter.com/MicrobiomDigest/status/1418944544725471233';...
    'https://twitter.com/MicrobiomDigest/status/1421128647554043916';...
    'https://twitter.com/MicrobiomDigest/status/1421906162979586052';...
    'https://twitter.com/MicrobiomDigest/status/1425509092752248840';...
    'https://twitter.com/MicrobiomDigest/status/1427349252007682049';...
    'https://twitter.com/MicrobiomDigest/status/1428440750971367425';...
    'https://twitter.com/MicrobiomDigest/status/1433868993509593088';...
    'https://twitter.com/MicrobiomDigest/status/1434573897819697155';...
    'https://twitter.com/MicrobiomDigest/status/1434592598786727938';...
    'https://twitter.com/MicrobiomDigest/status/1417044698582163457';...
    };

%% 
originalImage={'https://pbs.twimg.com/media/E-Tbo4sXMAgQt-0?format=jpg&name=large';...
    'https://pbs.twimg.com/media/E6cNeyaWUAMTR2b?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E6lh2ZUWUAcsHrh?format=png&name=small';...
    'https://pbs.twimg.com/media/E6lk-yDWYAYNxyc?format=png&name=900x900';...
    'https://pbs.twimg.com/media/E6uWODsWYBcbvrc?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E7EaWY-XEAIccd-?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E7jc77zXsAQJGkY?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E7ugG-eXIAU-mI9?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E8hs9P8XIAg97MH?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E872j3LXoAQpbJA?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E9LXRJfX0AA1E8W?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E-YgK7jXoBEdcNB?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E-ihWWbX0AsG9gP?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E-iyIybWYAgvXeW?format=jpg&name=900x900';...
    'https://pbs.twimg.com/media/E6panPYWEAAjOUY?format=jpg&name=900x900';...
    };



%% Read image
% read directly from the twitter feed

test_image=imread(originalImage{1});
% if this does not work download and read locally
% test_image=imread('Desktop\test.jpg');
figure(1)
imagesc(test_image)
% get dimensions of image
[rows,cols,levels]=size(test_image);

% Detect panels
% convert to gray level
test_image_gray = rgb2gray(test_image);
% panels are contained within white regions, those above ~240 in all RGB channels
panel_regions =imopen(((sum(test_image,3)/3)<242),ones(3));
%figure(2)
%imagesc(panel_regions)
%
% now remove the details of the labels to get individual panels
[panel_regions_L,numRegions] = bwlabel(panel_regions);
panel_regions_P              = regionprops(panel_regions_L,'Area');
small_large_threshold       = max([panel_regions_P.Area])*otsuthresh([panel_regions_P.Area]/max([panel_regions_P.Area]));
disp([mean([panel_regions_P.Area]) median([panel_regions_P.Area]) small_large_threshold])
[panel_regions_large,numLarge]    = bwlabel(ismember(panel_regions_L,find([panel_regions_P.Area]>small_large_threshold)));

%imagesc(panel_regions_large+panel_regions)

% close the panels separately so not to merge them
panel_regions_large2=zeros(rows,cols);
for k=1:numLarge
    panel_regions_large2=panel_regions_large2+imerode(imfill(imclose(panel_regions_large==k,ones(129)),'holes'),ones(16));
end

%imagesc(panel_regions_large2*2+panel_regions)
%
% the erosion is important to remove edges
%[panel_individual,numPanels] = bwlabel(imerode(panel_regions_large2,ones(16))) ;
[panel_individual,numPanels] = bwlabel(panel_regions_large2) ;
%imagesc(panel_individual+panel_regions)
%
% test image to be examined, use edges as this reveals better the overlaps
edges_panels = edge(test_image_gray,'canny',[],2).*(panel_individual>0);
figure(2)
imagesc(panel_individual.*(1-edges_panels))
% Get the bounding boxes
panel_boundingBox   = regionprops(panel_individual,'BoundingBox');
boxes               = reshape([panel_boundingBox.BoundingBox],4,numPanels)';
minBoxes            = min(boxes);
[y1,x1]             = hist(round(boxes(:,2)/5),20);
numColsPanels       = max(y1);
[y1,x1]             = hist(round(boxes(:,1)/5),20);
numRowsPanels       = max(y1);


% clear for previous runs
clear cross_corr_panels new_image
new_image(rows,cols,numPanels)=0;
% Obtain the cross correlation of all panels
% there may be a better way to run this but for the time it works
for counter1=1:numPanels
    coords_1 = ceil(panel_boundingBox(counter1).BoundingBox);
    rows_1   = coords_1(2)+5:min(rows,coords_1(2)+minBoxes(4)-5);
    cols_1   = coords_1(1)+5:min(cols,coords_1(1)+minBoxes(3)-5);
    for counter2=1:numPanels
        disp([counter1 counter2])
        coords_2 = ceil(panel_boundingBox(counter2).BoundingBox);
        rows_2   = coords_2(2)+5:coords_2(2)+minBoxes(4)-5;
        cols_2   = coords_2(1)+5:coords_2(1)+minBoxes(3)-5;
        cross_corr_panels(:,:,counter1,counter2) = xcorr2(edges_panels(rows_1,cols_1),...
            edges_panels(rows_2,cols_2));
        new_image(rows_2(1:end-1),cols_2(1:end-1),counter1) =...
            cross_corr_panels(1:2:end-1,1:2:end-1,counter1,counter2)+...;
            cross_corr_panels(2:2:end,2:2:end,counter1,counter2);
    end
end
% % individual matches can be visualised
% figure(11)
% subplot(131)
% imagesc(edges_panels(rows_1,cols_1));
% subplot(132)
% imagesc(edges_panels(rows_2,cols_2));
% subplot(133)
% mesh(new_image(rows_2,cols_2,counter1));colorbar
% colormap hot
%     axis ij
% 
%     axis tight
% display as a montage of all correlations

%pos1=reshape([1:28],4,7);
% pos2=reshape(1:numPanels,numCols,numRows)';
% 
% figure(9)
% for kk=1:numPanels
%     disp(kk)
%     h{kk}=subplot(4,7, pos2(kk));
%     %kk=21;
%     mesh((new_image(:,:,kk)/max(max(new_image(:,:,kk))))  );
%     axis ij
%     axis off
%     axis tight
% end
%     colormap hot
% Change dimensions for better visualisation
% for kk=1:numPanels
%     %h{kk}.Position(1)=h{kk}.Position(1)+0.02;
%     %h{kk}.Position(2)=h{kk}.Position(2)-0.02;
%     h{kk}.Position(3)=0.09;
%     h{kk}.Position(4)=0.19;
%     %h{kk}.Visible='on';
%     %h{kk}.Color='k';
% end
    
    
% Display results as plots that show the cross correlation between panels
figure(9)
clf
pos2=reshape(1:numPanels,numColsPanels,numRowsPanels)';
for kk=1:numPanels
    disp(kk)
    currentPanelData        = cross_corr_panels(:,:,kk,:);
    currentPanelPeaks       = squeeze(max(max(currentPanelData)));
    currentPanelMainPeak    = currentPanelPeaks(kk);
    currentPanelOtherPeaks  = currentPanelPeaks(setdiff(1:numPanels,kk));
    secondMaxRemoved        = currentPanelOtherPeaks(currentPanelOtherPeaks~=max(currentPanelOtherPeaks));
    meanOtherPeaks          = mean(secondMaxRemoved);
    stdOtherPeaks           = std(secondMaxRemoved);
    upperLimit              = meanOtherPeaks+3*stdOtherPeaks+1;
%     meanOtherPeaks          = mean(currentPanelOtherPeaks);
%     stdOtherPeaks           = std(currentPanelOtherPeaks);
%     upperLimit              = meanOtherPeaks+3*stdOtherPeaks;
    
    
    peaksAboveUpper         = find(currentPanelPeaks>upperLimit);
    
    peaksAboveUpper_L       = find(currentPanelPeaks>upperLimit);
    peaksAboveUpper_V       = currentPanelPeaks(peaksAboveUpper_L);
    
    h{kk}=subplot(numRowsPanels,numColsPanels, pos2(kk));
    %kk=21;
    
    plot(1:numPanels,currentPanelPeaks,'k',kk,currentPanelMainPeak,'b*',[1 numPanels],[upperLimit upperLimit],'r:',peaksAboveUpper_L,peaksAboveUpper_V,'dm')
    %axis ij
    %axis off
    if numel(peaksAboveUpper_L)>1
        title(num2str(peaksAboveUpper_L),'fontsize',11)
    else
        title(num2str(peaksAboveUpper_L),'fontsize',9)
    end
    axis tight
    set(gca,'yscale','linear')
    grid on
end


%% Change dimensions for better visualisation
for kk=1:numPanels
    %h{kk}.Position(1)=h{kk}.Position(1)+0.02;
    %h{kk}.Position(2)=h{kk}.Position(2)-0.02;
    h{kk}.Position(3)=0.09;
    h{kk}.Position(4)=0.19;
    %h{kk}.Visible='on';
    %h{kk}.Color='k';
end


%% Display separately
h11=figure(11)
kk=21;
    mesh((new_image(:,:,kk)/max(max(new_image(:,:,kk))))  );
    axis ij
    axis tight
     set(gca,'Color','k')
     colormap hot
     
     %%
