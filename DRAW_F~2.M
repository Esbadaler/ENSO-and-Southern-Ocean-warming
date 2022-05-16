clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Calculate ENSO-teleconnect tauu and skewness of tauu
if exist('Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% 1.1 Calculate ENSO-teleconnect change in tauu
    load /datastore/wan143/Reanalysis_DATA/HadISST/Updated/Cal_Nino12_Nino34_Nino3_Nino4_monthly_HadiSST_187001_201912.mat
    sst_nino34=sst_nino34((1947-1870+1)*12+1:end);%%%%%%%%194801-201912
    [d_trend,d_detrend]=quadratic_trend(sst_nino34);
    sst_nino34=d_detrend./std(d_detrend);
    sst_nino34_06=sst_nino34(6:12:end);
    sst_nino34_06=sst_nino34_06(1:end-1); %%%%%%1948-2018
    sst_nino34_07=sst_nino34(7:12:end);
    sst_nino34_07=sst_nino34_07(1:end-1); %%%%%%1948-2018
    sst_nino34_08=sst_nino34(8:12:end);
    sst_nino34_08=sst_nino34_08(1:end-1); %%%%%%1948-2018
    sst_nino34_09=sst_nino34(9:12:end);
    sst_nino34_09=sst_nino34_09(1:end-1); %%%%%%1948-2018
    sst_nino34_10=sst_nino34(10:12:end);
    sst_nino34_10=sst_nino34_10(1:end-1); %%%%%%1948-2018
    sst_nino34_11=sst_nino34(11:12:end);
    sst_nino34_11=sst_nino34_11(1:end-1); %%%%%%1948-2018
    sst_nino34_12=sst_nino34(12:12:end);
    sst_nino34_12=sst_nino34_12(1:end-1); %%%%%%1948-2018
    sst_nino34_01=sst_nino34(1:12:end);
    sst_nino34_01=sst_nino34_01(2:end); %%%%%%1949-2019
    sst_nino34_02=sst_nino34(2:12:end);
    sst_nino34_02=sst_nino34_02(2:end); %%%%%%1949-2019
    sst_nino34_03=sst_nino34(3:12:end);
    sst_nino34_03=sst_nino34_03(2:end); %%%%%%1949-2019
    sst_nino34_04=sst_nino34(4:12:end);
    sst_nino34_04=sst_nino34_04(2:end); %%%%%%1949-2019
    sst_nino34_05=sst_nino34(5:12:end);
    sst_nino34_05=sst_nino34_05(2:end); %%%%%%1949-2019
    sst_nino34_DJF  =(sst_nino34_12+sst_nino34_01+sst_nino34_02)./3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sst_nino34_DJF  =sst_nino34_DJF(1948-1948+1:end); %%%%%%1948-2018
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file='/datastore/wan143/Reanalysis_DATA/NCEP/Updated_wind_stress/tauu_monthly_194801_201912_r360x180.mat';
    load(file)
    tauu=tauu((1947-1948+1)*12+1:(2019-1948+1)*12,:,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sst_mask=nc_varget('/datastore/wan143/Reanalysis_DATA/HadISST/Updated/Updated_To2019Dec/HadISST_sst_187001_201912_r360x180.nc','sst');
    sst_mask=squeeze(mean(sst_mask,1));
    for tt=1:size(tauu,1)
        tauu_d=squeeze(tauu(tt,:,:));
        tauu_d(find(isnan(sst_mask)==1))=nan;
        tauu(tt,:,:)=tauu_d;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_monthly_ano =nan(size(tauu));
    tauu_monthly_mean=nan(12,numel(lat),numel(lon));
    for mm=1:12
        tauu_monthly_mean(mm,:,:)=squeeze(mean(tauu(mm:12:end,:,:),1));
    end
    for tt=1:size(tauu,1)
        mm=mod(tt,12);
        if mm~=0
            tauu_monthly_ano(tt,:,:)=squeeze(tauu(tt,:,:))-squeeze(tauu_monthly_mean(mm,:,:));
        else
            tauu_monthly_ano(tt,:,:)=squeeze(tauu(tt,:,:))-squeeze(tauu_monthly_mean(12,:,:));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_monthly_detrend=nan(size(tauu_monthly_ano));
    for latt=1:numel(lat)
        for lonn=1:numel(lon)
            d=squeeze(tauu_monthly_ano(:,latt,lonn));
            if numel(find(isnan(d)))>0
            else
                [d_trend,d_detrend]=quadratic_trend(d);
                tauu_monthly_detrend(:,latt,lonn)=d_detrend;
            end
        end
    end
    tauu=tauu_monthly_detrend;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear tauu_monthly_detrend; clear tauu_monthly_ano
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tauu_J0M1_onto_nino34_DJF   =nan(numel(lat),numel(lon));
    cor_tauu_J0M1_onto_nino34_DJF     =nan(numel(lat),numel(lon));
    p_tauu_J0M1_onto_nino34_DJF       =nan(numel(lat),numel(lon));
    tauu_J0M1_EN                      =nan(numel(lat),numel(lon));
    tauu_J0M1_LN                      =nan(numel(lat),numel(lon));
    tauu_DJF_EN                       =nan(numel(lat),numel(lon));
    tauu_DJF_LN                       =nan(numel(lat),numel(lon));
    for latt=1:numel(lat)
        for lonn=1:numel(lon)
            d=squeeze(tauu(:,latt,lonn));
            if numel(find(isnan(d)))>0
            else
                d_06=d(6:12:end);
                d_06=d_06(1:end-1); %%%%%%1979-2018
                d_07=d(7:12:end);
                d_07=d_07(1:end-1); %%%%%%1979-2018
                d_08=d(8:12:end);
                d_08=d_08(1:end-1); %%%%%%1979-2018
                d_09=d(9:12:end);
                d_09=d_09(1:end-1); %%%%%%1979-2018
                d_10=d(10:12:end);
                d_10=d_10(1:end-1); %%%%%%1979-2018
                d_11=d(11:12:end);
                d_11=d_11(1:end-1); %%%%%%1979-2018
                d_12=d(12:12:end);
                d_12=d_12(1:end-1); %%%%%%1979-2018
                d_01=d(1:12:end);
                d_01=d_01(2:end); %%%%%%1980-2019
                d_02=d(2:12:end);
                d_02=d_02(2:end); %%%%%%1980-2019
                d_03=d(3:12:end);
                d_03=d_03(2:end); %%%%%%1980-2019
                d_04=d(4:12:end);
                d_04=d_04(2:end); %%%%%%1980-2019
                d_05=d(5:12:end);
                d_05=d_05(2:end); %%%%%%1980-2019
                d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
                d_DJF         =(d_12+d_01+d_02)./3;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sliding_window=11;
                half_window   =(sliding_window-1)/2;
                %%%%%%%%%%
                d_JJASONDJFMAM_ave=nan(size(d_JJASONDJFMAM));
                for yy=half_window+1:numel(d_JJASONDJFMAM)-half_window
                    d_JJASONDJFMAM_ave(yy)=mean(d_JJASONDJFMAM(yy-half_window:yy+half_window));
                end
                d_JJASONDJFMAM=d_JJASONDJFMAM-d_JJASONDJFMAM_ave;
                %%%%%%%%%%
                d_DJF_ave=nan(size(d_DJF));
                for yy=half_window+1:numel(d_DJF)-half_window
                    d_DJF_ave(yy)=mean(d_DJF(yy-half_window:yy+half_window));
                end
                d_DJF=d_DJF-d_DJF_ave;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [r rms slope interc n] = ew_regress(sst_nino34_DJF,d_JJASONDJFMAM);
                slope_tauu_J0M1_onto_nino34_DJF(latt,lonn)=slope;
                x=sst_nino34_DJF;
                y=d_JJASONDJFMAM;
                x(find(isnan(y)))=[];
                y(find(isnan(y)))=[];
                [r p]=my_corr(x(:),y(:));
                cor_tauu_J0M1_onto_nino34_DJF(latt,lonn)=r;
                p_tauu_J0M1_onto_nino34_DJF(latt,lonn)=p;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dd=d_JJASONDJFMAM(find(sst_nino34_DJF>=0.5));
                dd(find(isnan(dd)))=[];
                tauu_J0M1_EN(latt,lonn)=sum(dd);
                
                dd=d_JJASONDJFMAM(find(sst_nino34_DJF<=-0.5));
                dd(find(isnan(dd)))=[];
                tauu_J0M1_LN(latt,lonn)=sum(dd);
                
                dd=d_DJF(find(sst_nino34_DJF>=0.5));
                dd(find(isnan(dd)))=[];
                tauu_DJF_EN(latt,lonn)=sum(dd);
                
                dd=d_DJF(find(sst_nino34_DJF<=-0.5));
                dd(find(isnan(dd)))=[];
                tauu_DJF_LN(latt,lonn)=sum(dd);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification.mat lon lat slope_tauu_J0M1_onto_nino34_DJF cor_tauu_J0M1_onto_nino34_DJF p_tauu_J0M1_onto_nino34_DJF ...
         tauu_J0M1_EN tauu_J0M1_LN tauu_DJF_EN tauu_DJF_LN
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% 1.2 Calculate ENSO-teleconnect change in tauu, zonal ave
    tauu_1=nan(size(tauu,1),size(tauu,2));
    for tt=1:size(tauu,1)
        for latt=1:numel(lat)
            d=squeeze(tauu(tt,latt,:));
            if numel(find(isnan(d)))>=numel(d)*0.5
%             if numel(find(isnan(d)))==numel(d)
            else
                tauu_1(tt,latt)=mean_withoutnan(d);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tauu_J0M1_onto_nino34_DJF_1=nan(numel(lat),1);
    cor_tauu_J0M1_onto_nino34_DJF_1  =nan(numel(lat),1);
    p_tauu_J0M1_onto_nino34_DJF_1    =nan(numel(lat),1);
    %%%%%%%%%%%%%%%%
    year_nan=1948:2018;
    sliding_window=11;
    half_window   =(sliding_window-1)/2;
    for yy=1:numel(year_nan)
        if yy<half_window+1 | yy>numel(d_JJASONDJFMAM)-half_window
            year_nan(yy)=nan;
        end
    end 
    %%%%%%%%%%%%%%%%
    year_target=year_nan(find(sst_nino34_DJF>=0.5));
    year_target(find(isnan(year_target)))=[];
    tauu_J0M1_EN_1                   =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF<=-0.5));
    year_target(find(isnan(year_target)))=[];
    tauu_J0M1_LN_1                   =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF>=0.5));
    year_target(find(isnan(year_target)))=[];
    tauu_DJF_EN_1                    =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF<=-0.5));
    year_target(find(isnan(year_target)))=[];
    tauu_DJF_LN_1                    =nan(numel(year_target),numel(lat));
    
    
    year_target=year_nan(find(sst_nino34_DJF>=0.75));
    year_target(find(isnan(year_target)))=[];
    tauu_J0M1_EN0p75_1               =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF<=-0.75));
    year_target(find(isnan(year_target)))=[];
    tauu_J0M1_LN0p75_1               =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF>=0.75));
    year_target(find(isnan(year_target)))=[];
    tauu_DJF_EN0p75_1                =nan(numel(year_target),numel(lat));
    
    year_target=year_nan(find(sst_nino34_DJF<=-0.75));
    year_target(find(isnan(year_target)))=[];
    tauu_DJF_LN0p75_1                =nan(numel(year_target),numel(lat));
    
    for latt=1:numel(lat)
        d=squeeze(tauu_1(:,latt));
        if numel(find(isnan(d)))>0
        else
            d_06=d(6:12:end);
            d_06=d_06(1:end-1); %%%%%%1979-2018
            d_07=d(7:12:end);
            d_07=d_07(1:end-1); %%%%%%1979-2018
            d_08=d(8:12:end);
            d_08=d_08(1:end-1); %%%%%%1979-2018
            d_09=d(9:12:end);
            d_09=d_09(1:end-1); %%%%%%1979-2018
            d_10=d(10:12:end);
            d_10=d_10(1:end-1); %%%%%%1979-2018
            d_11=d(11:12:end);
            d_11=d_11(1:end-1); %%%%%%1979-2018
            d_12=d(12:12:end);
            d_12=d_12(1:end-1); %%%%%%1979-2018
            d_01=d(1:12:end);
            d_01=d_01(2:end); %%%%%%1980-2019
            d_02=d(2:12:end);
            d_02=d_02(2:end); %%%%%%1980-2019
            d_03=d(3:12:end);
            d_03=d_03(2:end); %%%%%%1980-2019
            d_04=d(4:12:end);
            d_04=d_04(2:end); %%%%%%1980-2019
            d_05=d(5:12:end);
            d_05=d_05(2:end); %%%%%%1980-2019
            d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
            d_DJF         =(d_12+d_01+d_02)./3;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sliding_window=11;
            half_window   =(sliding_window-1)/2;
            %%%%%%%%%%
            d_JJASONDJFMAM_ave=nan(size(d_JJASONDJFMAM));
            for yy=half_window+1:numel(d_JJASONDJFMAM)-half_window
                d_JJASONDJFMAM_ave(yy)=mean(d_JJASONDJFMAM(yy-half_window:yy+half_window));
            end
            d_JJASONDJFMAM=d_JJASONDJFMAM-d_JJASONDJFMAM_ave;
            %%%%%%%%%%
            d_DJF_ave=nan(size(d_DJF));
            for yy=half_window+1:numel(d_DJF)-half_window
                d_DJF_ave(yy)=mean(d_DJF(yy-half_window:yy+half_window));
            end
            d_DJF=d_DJF-d_DJF_ave;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [r rms slope interc n] = ew_regress(sst_nino34_DJF,d_JJASONDJFMAM);
            slope_tauu_J0M1_onto_nino34_DJF_1(latt)=slope;
            x=sst_nino34_DJF;
            y=d_JJASONDJFMAM;
            x(find(isnan(y)))=[];
            y(find(isnan(y)))=[];
            [r p]=my_corr(x(:),y(:));
            cor_tauu_J0M1_onto_nino34_DJF_1(latt)=r;
            p_tauu_J0M1_onto_nino34_DJF_1(latt)=p;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dd=d_JJASONDJFMAM(find(sst_nino34_DJF>=0.5));
            dd(find(isnan(dd)))=[];
            tauu_J0M1_EN_1(:,latt)=dd;
            
            dd=d_JJASONDJFMAM(find(sst_nino34_DJF<=-0.5));
            dd(find(isnan(dd)))=[];
            tauu_J0M1_LN_1(:,latt)=dd;
            
            dd=d_DJF(find(sst_nino34_DJF>=0.5));
            dd(find(isnan(dd)))=[];
            tauu_DJF_EN_1(:,latt)=dd;
            
            dd=d_DJF(find(sst_nino34_DJF<=-0.5));
            dd(find(isnan(dd)))=[];
            tauu_DJF_LN_1(:,latt)=dd;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dd=d_JJASONDJFMAM(find(sst_nino34_DJF>=0.75));
            dd(find(isnan(dd)))=[];
            tauu_J0M1_EN0p75_1(:,latt)=dd;
            
            dd=d_JJASONDJFMAM(find(sst_nino34_DJF<=-0.75));
            dd(find(isnan(dd)))=[];
            tauu_J0M1_LN0p75_1(:,latt)=dd;
            
            dd=d_DJF(find(sst_nino34_DJF>=0.75));
            dd(find(isnan(dd)))=[];
            tauu_DJF_EN0p75_1(:,latt)=dd;
            
            dd=d_DJF(find(sst_nino34_DJF<=-0.75));
            dd(find(isnan(dd)))=[];
            tauu_DJF_LN0p75_1(:,latt)=dd;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification.mat slope_tauu_J0M1_onto_nino34_DJF_1 cor_tauu_J0M1_onto_nino34_DJF_1 p_tauu_J0M1_onto_nino34_DJF_1 ...
         tauu_J0M1_EN_1 tauu_J0M1_LN_1 tauu_DJF_EN_1 tauu_DJF_LN_1 ...
         tauu_J0M1_EN0p75_1 tauu_J0M1_LN0p75_1 tauu_DJF_EN0p75_1 tauu_DJF_LN0p75_1 -append
    
else
    load Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification.mat 
end

if exist('Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification_tauu_clim.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file='/datastore/wan143/Reanalysis_DATA/NCEP/Updated_wind_stress/tauu_monthly_194801_201912_r360x180.mat';
    load(file)
    tauu=tauu((1947-1948+1)*12+1:(2019-1948+1)*12,:,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sst_mask=nc_varget('/datastore/wan143/Reanalysis_DATA/HadISST/Updated/Updated_To2019Dec/HadISST_sst_187001_201912_r360x180.nc','sst');
    sst_mask=squeeze(mean(sst_mask,1));
    for tt=1:size(tauu,1)
        tauu_d=squeeze(tauu(tt,:,:));
        tauu_d(find(isnan(sst_mask)==1))=nan;
        tauu(tt,:,:)=tauu_d;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_J0M1_clim=nan(numel(lat),numel(lon));
    for latt=1:numel(lat)
        for lonn=1:numel(lon)
            d=squeeze(tauu(:,latt,lonn));
            if numel(find(isnan(d)))>0
            else
                d_06=d(6:12:end);
                d_06=d_06(1:end-1); %%%%%%1979-2018
                d_07=d(7:12:end);
                d_07=d_07(1:end-1); %%%%%%1979-2018
                d_08=d(8:12:end);
                d_08=d_08(1:end-1); %%%%%%1979-2018
                d_09=d(9:12:end);
                d_09=d_09(1:end-1); %%%%%%1979-2018
                d_10=d(10:12:end);
                d_10=d_10(1:end-1); %%%%%%1979-2018
                d_11=d(11:12:end);
                d_11=d_11(1:end-1); %%%%%%1979-2018
                d_12=d(12:12:end);
                d_12=d_12(1:end-1); %%%%%%1979-2018
                d_01=d(1:12:end);
                d_01=d_01(2:end); %%%%%%1980-2019
                d_02=d(2:12:end);
                d_02=d_02(2:end); %%%%%%1980-2019
                d_03=d(3:12:end);
                d_03=d_03(2:end); %%%%%%1980-2019
                d_04=d(4:12:end);
                d_04=d_04(2:end); %%%%%%1980-2019
                d_05=d(5:12:end);
                d_05=d_05(2:end); %%%%%%1980-2019
                d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
                tauu_J0M1_clim(latt,lonn)=mean_withoutnan(d_JJASONDJFMAM);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification_tauu_clim.mat lon lat tauu_J0M1_clim
else
    load Data_for_Draw_figure_2_NCEP_Final_v2_1948_2018_11yr_for_rectification_tauu_clim.mat
end
    
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Draw test figure
figure
load mmjet.mat

subplot(3,2,1)
lon_left=0;
lon_right=360;
lat_lower=-70;
lat_upper=-20;

lon_1=lon(find(lon>=lon_left & lon<=lon_right));
lat_1=lat(find(lat>=lat_lower & lat<=lat_upper));
slope_tauu_J0M1_onto_nino34_DJF_d=slope_tauu_J0M1_onto_nino34_DJF(find(lat>=lat_lower & lat<=lat_upper),find(lon>=lon_left & lon<=lon_right));
[c h]=contourf(lon_1,lat_1,slope_tauu_J0M1_onto_nino34_DJF_d,[-0.1:0.001:0.1]);
set(h,'color','none')
caxis([-0.021 0.021])
colorbar('southoutside','Ticks',[-0.018 -0.012 -0.006 0 0.006 0.012 0.018])
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
ylabel('Latitude')
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
set(gca,'xtick',[0:30:360],'xticklabel',{'','','60^oE','','','','180^oW','','','','60^oW','',''})
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN','80^oN'})
% set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','','','40^oS','','','','0^o','','','','40^oN','','','','80^oN'})
text(lon_left,lat_upper+5,'{\bfA} ENSO induced zonal wind anomalies','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
set(gca,'position',[0.1300    0.7718    0.2512    0.1870])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_tauu_J0M1_onto_nino34_DJF_d=p_tauu_J0M1_onto_nino34_DJF(find(lat>=lat_lower & lat<=lat_upper),find(lon>=lon_left & lon<=lon_right));
for latt=1:3:numel(lat_1)
    for lonn=1:3:numel(lon_1)
        if p_tauu_J0M1_onto_nino34_DJF_d(latt,lonn)<0.1
            plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tauu_J0M1_clim_d=tauu_J0M1_clim(find(lat>=lat_lower & lat<=lat_upper),find(lon>=lon_left & lon<=lon_right));
% contour(lon_1,lat_1,tauu_J0M1_clim_d,[-0.3:0.05:-0.1],'color','b','linewidth',1)
% contour(lon_1,lat_1,tauu_J0M1_clim_d,[0.1:0.05:0.3],'color','r','linewidth',1)
% contour(lon_1,lat_1,tauu_J0M1_clim_d,[0 0],'color','k','linewidth',1)
% 


subplot(3,2,2)
slope_tauu_J0M1_onto_nino34_DJF_1_d=slope_tauu_J0M1_onto_nino34_DJF_1(find(lat>=lat_lower & lat<=lat_upper));
plot(slope_tauu_J0M1_onto_nino34_DJF_1_d,lat_1,'color','k','linestyle','-','linewidth',1)
set(gca,'tickdir','out')
set(gca,'fontsize',15)
hold on
text(-0.004,-15,'{\bfB} Zonal average','fontsize',15)
xlim([-0.004 0.004])
set(gca,'xtick',[-0.003:0.003:0.003],'xticklabel',[-0.003:0.003:0.003])
ylim([lat_lower lat_upper])
% set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN','80^oN'})
set(gca,'ytick',[-80:10:80],'yticklabel',{'','','','','','','','','','','','','','','','',''})
% set(gca,'xminortick','on','yminortick','on')
plot(zeros(numel(lat_1),1),lat_1,'k--')
set(gca,'position',[0.3917    0.7718    0.0703    0.1871])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_tauu_J0M1_onto_nino34_DJF_1_d=p_tauu_J0M1_onto_nino34_DJF_1(find(lat>=lat_lower & lat<=lat_upper));
% slope_tauu_J0M1_onto_nino34_DJF_1_d(find(p_tauu_J0M1_onto_nino34_DJF_1_d>0.1))=nan;
% plot(slope_tauu_J0M1_onto_nino34_DJF_1_d,lat_1,'color','k','linestyle','-','linewidth',3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tauu_J0M1_clim_1=nan(numel(lat),1);
% for latt=1:numel(lat)
%     d=squeeze(tauu_J0M1_clim(latt,:));
%     if numel(find(isnan(d)))>=numel(d)*0.5
%         %             if numel(find(isnan(d)))==numel(d)
%     else
%         tauu_J0M1_clim_1(latt)=mean_withoutnan(d);
%     end
% end
% tauu_J0M1_clim_1_d=tauu_J0M1_clim_1(find(lat>=lat_lower & lat<=lat_upper));
% plot(tauu_J0M1_clim_1_d,lat_1,'color','k','linestyle','-','linewidth',1)
% 




subplot(3,2,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_J0M1_EN0p75_1_d=tauu_J0M1_EN0p75_1(:,find(lat>=lat_lower & lat<=lat_upper));
tauu_J0M1_LN0p75_1_d=tauu_J0M1_LN0p75_1(:,find(lat>=lat_lower & lat<=lat_upper));
tauu_J0M1_recti0p75_1_d=sum(tauu_J0M1_EN0p75_1_d,1)+sum(tauu_J0M1_LN0p75_1_d,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_J0M1_EN_1_d=tauu_J0M1_EN_1(:,find(lat>=lat_lower & lat<=lat_upper));
tauu_J0M1_LN_1_d=tauu_J0M1_LN_1(:,find(lat>=lat_lower & lat<=lat_upper));
tauu_J0M1_recti_1_d=sum(tauu_J0M1_EN_1_d,1)+sum(tauu_J0M1_LN_1_d,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(lat_1,sum(tauu_J0M1_EN0p75_1_d,1),'color','r','linestyle','-','linewidth',2)
hold on
plot(lat_1,sum(tauu_J0M1_LN0p75_1_d,1),'color','b','linestyle','-','linewidth',1)
plot(lat_1,tauu_J0M1_recti0p75_1_d,'color','k','linestyle','-','linewidth',3)
% plot(lat_1,tauu_J0M1_recti_1_d,'color','k','linestyle','--','linewidth',3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'tickdir','out')
set(gca,'fontsize',15)
ylabel('N m^-^2')
box on
text(lat_lower,0.14,'{\bfD} Cumulative ENSO zonal-mean zonal wind anomalies','fontsize',15)
xlim([lat_lower lat_upper])
set(gca,'xtick',[-80:10:80],'xticklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN','80^oN'})
ylim([-0.12 0.12])
set(gca,'ytick',[-0.1:0.04:0.1],'yticklabel',[-0.1:0.04:0.1])
set(gca,'xminortick','on','yminortick','on')
plot(lat_1,zeros(numel(lat_1),1),'k--','linewidth',3)
legend('Cumulative El Nino','Cumulative La Nina','Residual')
% legend('Cumulative EN','Cumulative LN','Residual','Residual(0.5)')
set(gca,'position',[0.1300    0.0453    0.3346    0.2200])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.1 Load in stoke streamfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Climatological SF and ENSO-SF
% file='/datastore/wan143/Reanalysis_DATA/NCEP/Calculate_stokes_streamfunction_from_vwnd_194801_201912_r360x180.mat';
load Calculate_stokes_streamfunction_from_vwnd_194801_201912_r360x180_AccordingToSurfacePressure_v6.mat
stoke_sf=stoke_sf((1947-1948+1)*12+1:end,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.2 Load in Nino3.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Cal_Nino12_Nino34_Nino3_Nino4_monthly_HadiSST_187001_201912.mat
sst_nino34=sst_nino34((1947-1870+1)*12+1:end);%%%%%%%%194801-201912
[d_trend,d_detrend]=quadratic_trend(sst_nino34);
sst_nino34=d_detrend./std(d_detrend);
sst_nino34_06=sst_nino34(6:12:end);
sst_nino34_06=sst_nino34_06(1:end-1); %%%%%%1948-2018
sst_nino34_07=sst_nino34(7:12:end);
sst_nino34_07=sst_nino34_07(1:end-1); %%%%%%1948-2018
sst_nino34_08=sst_nino34(8:12:end);
sst_nino34_08=sst_nino34_08(1:end-1); %%%%%%1948-2018
sst_nino34_09=sst_nino34(9:12:end);
sst_nino34_09=sst_nino34_09(1:end-1); %%%%%%1948-2018
sst_nino34_10=sst_nino34(10:12:end);
sst_nino34_10=sst_nino34_10(1:end-1); %%%%%%1948-2018
sst_nino34_11=sst_nino34(11:12:end);
sst_nino34_11=sst_nino34_11(1:end-1); %%%%%%1948-2018
sst_nino34_12=sst_nino34(12:12:end);
sst_nino34_12=sst_nino34_12(1:end-1); %%%%%%1948-2018
sst_nino34_01=sst_nino34(1:12:end);
sst_nino34_01=sst_nino34_01(2:end); %%%%%%1949-2019
sst_nino34_02=sst_nino34(2:12:end);
sst_nino34_02=sst_nino34_02(2:end); %%%%%%1949-2019
sst_nino34_03=sst_nino34(3:12:end);
sst_nino34_03=sst_nino34_03(2:end); %%%%%%1949-2019
sst_nino34_04=sst_nino34(4:12:end);
sst_nino34_04=sst_nino34_04(2:end); %%%%%%1949-2019
sst_nino34_05=sst_nino34(5:12:end);
sst_nino34_05=sst_nino34_05(2:end); %%%%%%1949-2019
sst_nino34_DJF  =(sst_nino34_12+sst_nino34_01+sst_nino34_02)./3;
sst_nino34_ONDJF=(sst_nino34_10+sst_nino34_11+sst_nino34_12+sst_nino34_01+sst_nino34_02)./5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.4 Regress SF onto Nino3.4
stoke_sf_ano =nan(size(stoke_sf));
stoke_sf_mean=nan(12,numel(lev_sf),numel(lat));
for mm=1:12
    stoke_sf_mean(mm,:,:)=squeeze(mean(stoke_sf(mm:12:end,:,:),1));
end
for tt=1:size(stoke_sf,1)
    mm=mod(tt,12);
    if mm~=0
        stoke_sf_ano(tt,:,:)=squeeze(stoke_sf(tt,:,:))-squeeze(stoke_sf_mean(mm,:,:));
    else
        stoke_sf_ano(tt,:,:)=squeeze(stoke_sf(tt,:,:))-squeeze(stoke_sf_mean(12,:,:));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stoke_sf_detrend=nan(size(stoke_sf_ano));
for levv=1:numel(lev_sf)
    for latt=1:numel(lat)
        d=squeeze(stoke_sf_ano(:,levv,latt));
        if numel(find(isnan(d)))>0
        else
            [d_trend,d_detrend]=quadratic_trend(d);
            stoke_sf_detrend(:,levv,latt)=d_detrend;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_StokeSF_J0M1_onto_nino34_DJF   =nan(numel(lev_sf),numel(lat));
slope_StokeSF_J0M1_onto_nino34_ONDJF =nan(numel(lev_sf),numel(lat));
cor_StokeSF_J0M1_onto_nino34_DJF     =nan(numel(lev_sf),numel(lat));
cor_StokeSF_J0M1_onto_nino34_ONDJF   =nan(numel(lev_sf),numel(lat));
p_StokeSF_J0M1_onto_nino34_DJF       =nan(numel(lev_sf),numel(lat));
p_StokeSF_J0M1_onto_nino34_ONDJF     =nan(numel(lev_sf),numel(lat));
StokeSF_J0M1_EN_0p5                  =nan(numel(lev_sf),numel(lat));
StokeSF_J0M1_LN_0p5                  =nan(numel(lev_sf),numel(lat));
StokeSF_J0M1_EN_0p75                 =nan(numel(lev_sf),numel(lat));
StokeSF_J0M1_LN_0p75                 =nan(numel(lev_sf),numel(lat));
for levv=1:numel(lev_sf)
    for latt=1:numel(lat)
        d=squeeze(stoke_sf_detrend(:,levv,latt));
        if numel(find(isnan(d)))>0
        else
            d_06=d(6:12:end);
            d_06=d_06(1:end-1); %%%%%%1979-2018
            d_07=d(7:12:end);
            d_07=d_07(1:end-1); %%%%%%1979-2018
            d_08=d(8:12:end);
            d_08=d_08(1:end-1); %%%%%%1979-2018
            d_09=d(9:12:end);
            d_09=d_09(1:end-1); %%%%%%1979-2018
            d_10=d(10:12:end);
            d_10=d_10(1:end-1); %%%%%%1979-2018
            d_11=d(11:12:end);
            d_11=d_11(1:end-1); %%%%%%1979-2018
            d_12=d(12:12:end);
            d_12=d_12(1:end-1); %%%%%%1979-2018
            d_01=d(1:12:end);
            d_01=d_01(2:end); %%%%%%1980-2019
            d_02=d(2:12:end);
            d_02=d_02(2:end); %%%%%%1980-2019
            d_03=d(3:12:end);
            d_03=d_03(2:end); %%%%%%1980-2019
            d_04=d(4:12:end);
            d_04=d_04(2:end); %%%%%%1980-2019
            d_05=d(5:12:end);
            d_05=d_05(2:end); %%%%%%1980-2019
            d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            [r rms slope interc n] = ew_regress(sst_nino34_DJF,d_JJASONDJFMAM);
            slope_StokeSF_J0M1_onto_nino34_DJF(levv,latt)=slope;
            x=sst_nino34_DJF;
            y=d_JJASONDJFMAM;
            [r p]=my_corr(x(:),y(:));
            cor_StokeSF_J0M1_onto_nino34_DJF(levv,latt)=r;
            p_StokeSF_J0M1_onto_nino34_DJF(levv,latt)=p;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            [r rms slope interc n] = ew_regress(sst_nino34_ONDJF,d_JJASONDJFMAM);
            slope_StokeSF_J0M1_onto_nino34_ONDJF(levv,latt)=slope;
            x=sst_nino34_ONDJF;
            y=d_JJASONDJFMAM;
            [r p]=my_corr(x(:),y(:));
            cor_StokeSF_J0M1_onto_nino34_ONDJF(levv,latt)=r;
            p_StokeSF_J0M1_onto_nino34_ONDJF(levv,latt)=p;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EN & LN
            StokeSF_J0M1_EN_0p5(levv,latt)=mean(d_JJASONDJFMAM(find(sst_nino34_DJF>=0.5)));
            StokeSF_J0M1_LN_0p5(levv,latt)=mean(d_JJASONDJFMAM(find(sst_nino34_DJF<=-0.5)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EN & LN
            StokeSF_J0M1_EN_0p75(levv,latt)=mean(d_JJASONDJFMAM(find(sst_nino34_DJF>=0.75)));
            StokeSF_J0M1_LN_0p75(levv,latt)=mean(d_JJASONDJFMAM(find(sst_nino34_DJF<=-0.75)));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.5 Draw regression pattern
load mmjet.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,3)

lat_lower=-80;
lat_upper=20;
lev_lower=0;
lev_upper=100000;

lat_1=lat(find(lat>=lat_lower & lat<=lat_upper));
lev_1=lev_sf(find(lev_sf>=lev_lower & lev_sf<=lev_upper));

% cor_StokeSF_J0M1_onto_nino34_DJF(find(p_StokeSF_J0M1_onto_nino34_DJF>0.1))=nan;
slope_StokeSF_J0M1_onto_nino34_DJF_d=slope_StokeSF_J0M1_onto_nino34_DJF(find(lev_sf>=lev_lower & lev_sf<=lev_upper),find(lat>=lat_lower & lat<=lat_upper));
[c h]=contourf(lat_1,lev_1,slope_StokeSF_J0M1_onto_nino34_DJF_d,[-14e9:0.1e9:14e9]);
hold on
set(h,'color','none')
caxis([-1.4 1.4].*1.0e9)
colorbar('southoutside','Ticks',[-1.2 -0.8 -0.4 0 0.4 0.8 1.2].*1.0e9)
colormap(gca,mm_jet_normal)
% caxis([-0.5 0.5])
% colorbar('Ticks',[-0.4 -0.2 0 0.2 0.4])
% colormap(gca,cmjet)
set(gca,'tickdir','out')
set(gca,'xtick',[-80:10:80],'xticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
set(gca,'ydir','reverse')
set(gca,'ytick',[10000:10000:100000],'yticklabel',{'100','','300','','500','','700','','900',''})
ylabel('mb')
text(lat_lower,lev_lower-10000,'{\bfC} ENSO-induced zonal mean circulation anomalies','fontsize',15)
hold on
xlim([lat_lower lat_upper])
ylim([lev_lower lev_upper])
set(gca,'fontsize',15)
set(gca,'xminortick','on','yminortick','on')
set(gca,'position',[0.1300    0.4234    0.3346    0.1870])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_StokeSF_J0M1_onto_nino34_DJF_d=p_StokeSF_J0M1_onto_nino34_DJF(find(lev_sf>=lev_lower & lev_sf<=lev_upper),find(lat>=lat_lower & lat<=lat_upper));
for levv=1:numel(lev_1)
    for latt=1:numel(lat_1)
        if p_StokeSF_J0M1_onto_nino34_DJF_d(levv,latt)<0.1
            plot(lat_1(latt),lev_1(levv),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',0.5)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stoke_sf_12=stoke_sf(12:12:end-1,:,:);
stoke_sf_13=stoke_sf(13:12:end-1,:,:);
stoke_sf_14=stoke_sf(14:12:end-1,:,:);
stoke_sf_djf=(stoke_sf_12+stoke_sf_13+stoke_sf_14)./3;
stoke_sf_djf_clim=squeeze(mean(stoke_sf_djf,1));
stoke_sf_djf_clim_d=stoke_sf_djf_clim(find(lev_sf>=lev_lower & lev_sf<=lev_upper),find(lat>=lat_lower & lat<=lat_upper));
% contour(lat_1,lev_1,stoke_sf_djf_clim_d,[-5e10 -4e10 -3e10 -2e10 -1e10 -0.5e10 -0.3e10 -0.1e10],'linestyle','--','color','k','linewidth',1);
% contour(lat_1,lev_1,stoke_sf_djf_clim_d,[17e10 13e10 9e10 5e10 4e10 3e10 2e10 1e10],'linestyle','-','color','k','linewidth',1);
contour(lat_1,lev_1,stoke_sf_djf_clim_d,[0 0],'linestyle','-','color','k','linewidth',3);
contour(lat_1,lev_1,stoke_sf_djf_clim_d,[-3e10 -1e10 -0.5e10 -0.1e10],'linestyle','--','color','k','linewidth',1);
contour(lat_1,lev_1,stoke_sf_djf_clim_d,[13e10 9e10 3e10 1e10 0.5e10 0.1e10],'linestyle','-','color','k','linewidth',1);



























