%% Read image
% read directly from the twitter feed
test_image=imread('https://pbs.twimg.com/media/E-Tbo4sXMAgQt-0?format=jpg&name=large');
% if this does not work download and read locally
% test_image=imread('Desktop\test.jpg');
imagesc(test_image)
% get dimensions of image
[rows,cols,levels]=size(test_image);

%% Detect panels
% convert to gray level
test_image_gray = rgb2gray(test_image);
% panels are contained within white regions, those above ~240 in all RGB channels
panel_regions =((sum(test_image,3)/3)<240);
imagesc(panel_regions)

% now remove the details of the labels to get individual panels
% the erosion is important to remove edges
[panel_individual,numPanels] = bwlabel(imerode(imopen(panel_regions,ones(30,30)),ones(16)) );
imagesc(panel_individual)

% test image to be examined, use edges as this reveals better the overlaps
test_panels = edge(test_image_gray,'canny',[],2).*(panel_individual>0);
imagesc(test_panels)
% Get the bounding boxes
panel_boundingBox = regionprops(panel_individual,'BoundingBox');
boxes             = reshape([panel_boundingBox.BoundingBox],4,numPanels)';
minBoxes          = min(boxes);
% clear for previous runs
clear cross_corr_panels new_image
%% Obtain the cross correlation of all panels
% there may be a better way to run this but for the time it works
for counter1=1:numPanels
    coords_1 = ceil(panel_boundingBox(counter1).BoundingBox);
    rows_1   = coords_1(2)+5:coords_1(2)+minBoxes(4)-4;
    cols_1   = coords_1(1)+5:coords_1(1)+minBoxes(3)-4;
    for counter2=1:numPanels
        disp([counter1 counter2])
        coords_2 = ceil(panel_boundingBox(counter2).BoundingBox);
        rows_2   = coords_2(2)+5:coords_2(2)+minBoxes(4)-4;
        cols_2   = coords_2(1)+5:coords_2(1)+minBoxes(3)-4;
        cross_corr_panels(:,:,counter1,counter2) = xcorr2(test_panels(rows_1,cols_1),...
            test_panels(rows_2,cols_2));
        new_image(rows_2(1:end-1),cols_2(1:end-1),counter1) =...
            cross_corr_panels(1:2:end-1,1:2:end-1,counter1,counter2)+...;
            cross_corr_panels(2:2:end,2:2:end,counter1,counter2);
    end
end
%% individual matches can be visualised
% figure(1)
% imagesc(test_panels(rows_1,cols_1));colorbar
% figure(2)
% imagesc(test_panels(rows_2,cols_2));colorbar
% figure(3)
% mesh(new_image(rows_2,cols_2,counter1));colorbar
%%
% for counter1=1:numPanels
%     coords_1 = ceil(panel_boundingBox(counter1).BoundingBox);
%     rows_1   = coords_1(2)+5:coords_1(2)+minBoxes(4)-4;
%     cols_1   = coords_1(1)+5:coords_1(1)+minBoxes(3)-4;
%     for counter2=1:numPanels
%         disp([counter1 counter2])
%         coords_2 = ceil(panel_boundingBox(counter2).BoundingBox);
%         rows_2   = coords_2(2)+5:coords_2(2)+minBoxes(4)-4;
%         cols_2   = coords_2(1)+5:coords_2(1)+minBoxes(3)-4;
%         new_image(rows_2,cols_2,counter1) =cross_corr_panels(1:2:end,1:2:end,counter1,counter2);
%     end
% end

%%

pos1=reshape([1:28],4,7);
pos2=reshape([1:28],7,4)';

figure(8)
for kk=1:numPanels
    disp(kk)
    h{kk}=subplot(4,7, pos2(kk));
    %kk=21;
    mesh((new_image(:,:,kk)/max(max(new_image(:,:,kk))))  );
    axis ij
    axis off
    axis tight
end
    colormap hot
    set(gca,'Color','k')
%%
for kk=1:numPanels
    %h{kk}.Position(1)=h{kk}.Position(1)+0.02;
    %h{kk}.Position(2)=h{kk}.Position(2)-0.02;
    
    h{kk}.Position(3)=0.1;
    h{kk}.Position(4)=0.22;
    h{kk}.Visible='on';
    h{kk}.Color='k';
    
end


%%
figure(1)
kk=10;
    mesh((new_image(:,:,kk)/max(max(new_image(:,:,kk))))  );
    axis ij
    axis tight
     set(gca,'Color','k')
     
     %%