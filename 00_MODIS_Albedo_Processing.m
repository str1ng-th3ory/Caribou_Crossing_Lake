clc
clear

lc_file = 'landcover/MCD12Q1.061_LC_Type1_doy2011001_aid0001.tif';

%%generate water mask
lake_img = imread(lc_file);

%convert to 
water_mask = single(lake_img);
water_mask(water_mask ~= 17) = nan;
water_mask(water_mask == 17) = 1;
 
%imagesc(water_mask); show the mask
[rows,cols] = size(water_mask);

% %%rearrange different files to sub folders
% for iyear = 2001:2016
%    files = dir(['alb\*doy',num2str(iyear),'*.tif']); 
%    mkdir(['E:\alb\',num2str(iyear)]);
%    iyear
%    for ifile = 1:length(files)
%       [a,b] = movefile([files(ifile,1).folder,'\',files(ifile,1).name],['E:\alb\',num2str(iyear)]); 
%    end
% end


%%climatology calculation

alb_count_all = single(zeros(rows,cols,365)); 
alb_all = single(zeros(rows,cols,365));

for iyear = 2001 : 2022
   iyear
   files = dir(['alb/',num2str(iyear),'/MCD43A3.061_Albedo_BSA_shortwave_doy*.tif']); 
   
   alb_year = single(zeros(rows,cols,365));
   alb_count = single(zeros(rows,cols,365));
   for ifile = 1:length(files)
   
        day_alb = single(imread([files(ifile,1).folder,'/',files(ifile,1).name]))/1000;
        iday = str2num(files(ifile,1).name(41:43));
        if iday>365
           continue 
        end
        day_alb(day_alb<=0 | day_alb>=1)=nan;
        alb_year(:,:,iday) = day_alb .* water_mask; % only keep water pixel
   end
   alb_year(isnan(alb_year)) = 0;
   alb_count = alb_year; alb_count(alb_count~=0) = 1;
   
   alb_all = alb_all + alb_year;

   alb_count_all = alb_count_all + alb_count;
end

alb_clim = alb_all./alb_count_all; %calculate mean
alb_clim(alb_clim<=0 | alb_clim>=1) = nan; %% double check if the downloaded albedo has been converted to 0~1

for iyear = 2017:2022 
   iyear %%target year
   
   files = dir(['alb/',num2str(iyear),'/MCD43A3.061_Albedo_BSA_shortwave_doy*.tif']);
   alb_year = single(nan(rows,cols,365));
   for ifile = 1:length(files)
       day_alb = single(imread([files(ifile,1).folder,'/',files(ifile,1).name]))/1000;
       iday = str2num(files(ifile,1).name(41:43));
       day_alb(day_alb<=0 | day_alb>=1)=nan;
       alb_year(:,:,iday) = day_alb .* water_mask; % only keep water pixel
   end
   
   files = dir(['alb/',num2str(iyear),'/MCD43A3.061_BRDF_Albedo_Band_Mandatory_Quality_shortwave_doy*.tif']); % get good retrieval
   qa_year = single(zeros(rows,cols,365));
   for ifile = 1:length(files)
       day_qa = single(imread([files(ifile,1).folder,'/',files(ifile,1).name]));
       iday = str2num(files(ifile,1).name(65:67));
       qa_year(:,:,iday) = day_qa; % only keep water pixel
   end
   
   alb_target = alb_year;
   alb_target(qa_year>=2) = nan; %quality check
   
   output_alb = alb_target; %%set the clear-sky albedos in the output matrix, nan values included
   
   %%fill nan values using KF
   for m = 1:rows
       m
       for n = 1:cols
           if isnan(water_mask(m,n))
               continue
           end
           
           pixel_alb = squeeze(alb_target(m,n,:));
           a=find(~isnan(pixel_alb)); % find the first and last valid day, and give up the nearby values due to the low retrieval accuracy
           start_day = 75;%a(1);
           end_day = 295;%a(end);
           buffer_days = 15; % changeable
           pixel_alb(1:start_day+buffer_days) = nan;
           pixel_alb(end_day-buffer_days+1:end) = nan;
           
           pixel_clim =  squeeze(alb_clim(m,n,:));
           pixel_clim = fillmissing(pixel_clim,'nearest');% temporal interpolation, just make sure the climatology is continous
           pixel_clim(1:start_day+buffer_days) = nan;
           pixel_clim(end_day-buffer_days+1:end) = nan;
           
           
           %%filling, paper: 2.2.3 of Jia, Aolin, et al. "Improved cloudy-sky snow albedo estimates using passive microwave and VIIRS data." ...
           %% ISPRS Journal of Photogrammetry and Remote Sensing 196 (2023): 340-355.
           dynamic_alb = pixel_clim;
           dT_model = squeeze(dynamic_alb -[NaN;dynamic_alb(1:end-1)]); %difference between the day and the previous day (day-1)
           mark = pixel_alb;mark(~isnan(mark)) = 1;
           y_kalman = nan(365,1);%result
           P_fil = 3;%climatology error %relative ratio
           R = 3;%遥感观测误差 精度0.05 常数
           Q=P_fil;%初始值
           x_fil=[dynamic_alb(start_day+buffer_days+1)];
           
           for ii = start_day+buffer_days+2:end_day-buffer_days %filtering from the second day
               adjust_adj_clim = nan;
               if isnan(mark(ii))
                   win_m = [m-2:m+2]; win_m(win_m <1 | win_m > rows) = [];
                   win_n = [n-2:n+2]; win_n(win_n <1 | win_n > cols) = [];
                   adj_alb = squeeze(alb_target(win_m,win_n,ii));
                   if length(adj_alb(~isnan(adj_alb)))>2 %%spatial adjacent clear-sky pixels at least 3
                       adj_clim = squeeze(alb_clim(win_m,win_n,ii));
                       diff = adj_alb - adj_clim; %baise between this year and climatology year
                       adjust_adj_clim = alb_clim(m,n,ii) + median(diff(:),'omitnan');
                       adjust_adj_clim(adjust_adj_clim<=0 | adjust_adj_clim>=1) = nan;
                   end
               end
               A = [];%simulation model 误差传播
               if ii==start_day+buffer_days+2
                   A=[1 dT_model(start_day+buffer_days+2)./(dynamic_alb(start_day+buffer_days+1))];
               else
                   A=[1 dT_model(ii)./(y_kalman(ii-1))];
               end
               if isnan(mark(ii,1))
                   if ~isnan(adjust_adj_clim)% if spatially adjacient pixels are available
                       y_kalman(ii,1) = adjust_adj_clim;
                       x_fil=[y_kalman(ii,1)];
                       continue
                   end
                   if ii==start_day+buffer_days+2
                       y_kalman(ii,1) = sum(dynamic_alb(start_day+buffer_days+1) *(A));
                   else
                       y_kalman(ii,1) = sum(y_kalman(ii-1,1)*[1 dT_model(ii)./y_kalman(ii-1,1)]);
                   end
                   if y_kalman(ii,1) <= 0 | y_kalman(ii,1)>=1
                       y_kalman(ii,1) = mean((y_kalman(ii-3:ii-1,1) - pixel_clim(ii-3:ii-1)),'omitnan') + mean((pixel_clim(ii-3:ii-1)),'omitnan');
                   end
                   x_fil=[y_kalman(ii,1)];
                   P_fil=A*P_fil*A'+Q;
               else
                   x_pre=sum(A*(x_fil));
                   P_pre=A*P_fil*A'+Q; %acumulated climatology error (weight)
                   K=P_pre./(P_pre+R);%kalman gain
                   x_fil=x_pre+K*(pixel_alb(ii)-x_pre);
                   P_fil=(1-K)*P_pre;
                   y_kalman(ii,1)=x_fil(1);
               end
               
           end
           % fill the nan values by cloudy-sky estimations
           pixel_alb(isnan(pixel_alb)) = y_kalman(isnan(pixel_alb)); % only use the cloudy-sky pixels, y_kalman is the filtered results
           output_alb(m,n,:) = pixel_alb;%final output
       end
   end
   
   close all

   %for i = 1:365
   %    h=imagesc(output_alb(:,:,i));
   %    colorbar
   %    set(h,'alphadata',~isnan(output_alb(:,:,i)));
   %    title(num2str(i,'%03d'))
   %    caxis([0,0.7])
   %    pause(0.3)
   %    colormap jet
       
   %end
   
   save(['output_',num2str(iyear),'.mat'],'output_alb');

   %%test
   %h1=figure
   %m=211; n=450;
   %len = find(~isnan(output_alb(m,n,:)));
   %plot(len,squeeze(output_alb(m,n,len)),'LineWidth',3);hold on
   %plot(len,squeeze(alb_clim(m,n,len)),'LineWidth',2);
   %plot(len,squeeze(alb_year(m,n,len)),'.','MarkerSize',15)
   %legend('Final Output','Climatology','Clear-sky Only');
   %set(gcf,'position',[5 150 1000 800]);
   
   %h1=figure
   %m=100; n=200;
   %len = find(~isnan(output_alb(m,n,:)));
   %plot(len,squeeze(output_alb(m,n,len)),'LineWidth',3);hold on
   %plot(len,squeeze(alb_clim(m,n,len)),'LineWidth',2);
   %plot(len,squeeze(alb_year(m,n,len)),'.','MarkerSize',15)
   %legend('Final Output','Climatology','Clear-sky Only');
   %set(gcf,'position',[5 150 1000 800]);
   
   %h1=figure
   %m=330; n=630;
   %len = find(~isnan(output_alb(m,n,:)));
   %plot(len,squeeze(output_alb(m,n,len)),'LineWidth',3);hold on
   %plot(len,squeeze(alb_clim(m,n,len)),'LineWidth',2);
   %plot(len,squeeze(alb_year(m,n,len)),'.','MarkerSize',15)
   %legend('Final Output','Climatology','Clear-sky Only');
   %set(gcf,'position',[5 150 1000 800]);


end



%%
