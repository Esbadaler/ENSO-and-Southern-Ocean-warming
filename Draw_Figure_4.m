clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.1 Load in Nino SST index
load Calculate_Nino3_Nino4_Nino34_DJF_std_skewness_CMIP6.mat
models_sst=models;
clear models

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.2 Load in GMT using ST
load Calculate_Global_mean_ST_190001_209912_CMIP6.mat
models_ST=models;
clear models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pick_up_ST=1:37;
pick_up_ST([1 2 19])=[];
models_ST =models_ST(pick_up_ST);
ts_GMT_monthly=ts_GMT_monthly(pick_up_ST,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts_GMT_ave_s1=nan(numel(models_ST),1);
ts_GMT_ave_s2=nan(numel(models_ST),1);
for num=1:numel(models_ST)
    ts_GMT_ave_s1(num)=mean_withoutnan(ts_GMT_monthly(num,1:1200));
    ts_GMT_ave_s2(num)=mean_withoutnan(ts_GMT_monthly(num,1201:end));
end
ts_GMT_diff=ts_GMT_ave_s2-ts_GMT_ave_s1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Calculate mean state in ocean temperature, all seasons.
if exist('Data_for_Draw_relationship_change_in_thetaoAveGlobal_monthly_vs_ECS_all_depth.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Ocean temperature
    file_base=['/datastore/wan143/thetao/cmip6/monthly_historical_ssp585/1900-2099/thetao_Omon_' cell2mat(models_sst(1)) '_monthly_190001_209912.mat'];
    load(file_base)
    lev_base=[10  20  30  40  50  60  70  80  90  100 ...
              110 120 130 140 150 160 170 180 190 200 ...
              210 220 250 280 300 340 400 500 600 700 800 900 1000 ...
              1300 1500 1800 2000 2300 2500 2800 3000 ...
              3300 3500 3800 4000 4300 4500 4800 5000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    lon_left =0;
    lon_right=360;
    lat_lower=-70;
    lat_upper=60;
    lat_base =lat(find(lat>=lat_lower & lat<=lat_upper));
    thetaoGlobal_monthly_mean_s1=nan(numel(models_sst),numel(lev_base),numel(lat_base));
    thetaoGlobal_monthly_mean_s2=nan(numel(models_sst),numel(lev_base),numel(lat_base));
    clear thetao_monthly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/thetao/cmip6/monthly_historical_ssp585/1900-2099/thetao_Omon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load(file)
            if size(thetao_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetao_monthly_1=thetao_monthly(:,:,find(lat>=lat_lower & lat<=lat_upper),find(lon>=lon_left & lon<=lon_right));
            thetao_monthly_Global=nan(size(thetao_monthly_1,1),size(thetao_monthly_1,2),size(thetao_monthly_1,3));
            for tt=1:size(thetao_monthly_1,1)
                for levv=1:size(thetao_monthly_1,2)
                    for latt=1:size(thetao_monthly_1,3)
                        d=squeeze(thetao_monthly_1(tt,levv,latt,:));
                        thetao_monthly_Global(tt,levv,latt)=mean_withoutnan(d(:));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetao_monthly_Global_regrided=nan(size(thetao_monthly_Global,1),numel(lev_base),numel(lat_base));
            for tt=1:size(thetao_monthly_Global,1)
                for latt=1:numel(lat_base)
                    d=squeeze(thetao_monthly_Global(tt,:,latt));
                    if numel(find(isnan(d)))==numel(d)
                    else
                        thetao_monthly_Global_regrided(tt,:,latt)=interp1(lev,d,lev_base);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for levv=1:numel(lev_base)
                for latt=1:numel(lat_base)
                    d=squeeze(thetao_monthly_Global_regrided(:,levv,latt));
                    if numel(find(isnan(d)))>=600
                    else
                        d_s1=d(1:1200);
                        d_s2=d(1201:end);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        thetaoGlobal_monthly_mean_s1(num,levv,latt)=mean_withoutnan(d_s1);
                        thetaoGlobal_monthly_mean_s2(num,levv,latt)=mean_withoutnan(d_s2);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save Data_for_Draw_relationship_change_in_thetaoAveGlobal_monthly_vs_ECS_all_depth.mat models_sst lev_base lat_base ...
         thetaoGlobal_monthly_mean_s1 thetaoGlobal_monthly_mean_s2 
     
else
    load Data_for_Draw_relationship_change_in_thetaoAveGlobal_monthly_vs_ECS_all_depth.mat
end

load Data_for_Draw_supplementary_figure_1_ttest.mat

lat_thetao=lat_base;
lev_thetao=lev_base;
clear lev_base; clear lat_base 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.1 Calculate mean state of Hur (zonal mean) and regression pattern of PD Hur onto Nino3.4 (interannual, DJF)
if exist('Data_for_Draw_relationship_Tauu_with_Nino34_Seasons_all_CMIP6_New7_corrected_New_New_New.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_base=['/datastore/wan143/TAUU/cmip6/monthly_historical_ssp585/1900-2099/tauu_Amon_' cell2mat(models_sst(1)) '_monthly_190001_209912.mat'];
    load(file_base)
    lon_base=lon;
    lat_base=lat;
    clear tauu_monthly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tauu_J0M1_and_nino34_O0F1_PD=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tauu_J0M1_and_nino34_O0F1_PD  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tauu_J0M1_and_nino34_O0F1_PD    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    slope_tauu_J0M1_and_nino34_DJF_PD=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tauu_J0M1_and_nino34_DJF_PD  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tauu_J0M1_and_nino34_DJF_PD    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tauu_J0M1_and_nino34_O0F1_CC=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tauu_J0M1_and_nino34_O0F1_CC  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tauu_J0M1_and_nino34_O0F1_CC    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    slope_tauu_J0M1_and_nino34_DJF_CC=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tauu_J0M1_and_nino34_DJF_CC  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tauu_J0M1_and_nino34_DJF_CC    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_monthly_mean_s1=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_monthly_mean_s2=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_djf_mean_s1=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_djf_mean_s2=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/TAUU/cmip6/monthly_historical_ssp585/1900-2099/tauu_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            load(file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(tauu_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    d=squeeze(tauu_monthly(:,latt,lonn));
                    tauu_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    tauu_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    d_12=d(12:12:end-1);
                    d_13=d(13:12:end-1);
                    d_14=d(14:12:end-1);
                    d_djf=(d_12+d_13+d_14)./3;
                    tauu_djf_mean_s1(num,latt,lonn)=mean_withoutnan(d_djf(1:100));
                    tauu_djf_mean_s2(num,latt,lonn)=mean_withoutnan(d_djf(101:end));
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/TAUU/cmip6/monthly_historical_ssp585/1900-2099/tauu_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            load(file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(tauu_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Calculate monthly anomalies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tauu_ano =nan(size(tauu_monthly));
            tauu_mean=nan(12,numel(lat),numel(lon));
            for mm=1:12
                tauu_mean(mm,:,:)=squeeze(mean(tauu_monthly(mm:12:1200,:,:),1));
            end
            for tt=1:size(tauu_monthly,1)
                mm=mod(tt,12);
                if mm~=0
                    tauu_ano(tt,:,:)=squeeze(tauu_monthly(tt,:,:))-squeeze(tauu_mean(mm,:,:));
                else
                    tauu_ano(tt,:,:)=squeeze(tauu_monthly(tt,:,:))-squeeze(tauu_mean(12,:,:));
                end
            end
            tauu_monthly=tauu_ano;
            clear tauu_ano
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%Quadratic detrend
            tauu_monthly_detrend=nan(size(tauu_monthly));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tauu_monthly(:,latt,lonn));
                    if numel(find(isnan(d)))>=100
                    elseif numel(find(isnan(d)))>0 & numel(find(isnan(d)))<100
                        d_1=interp1_nan(d);
                        if numel(find(isnan(d_1)))>0
                            year=1:2400;
                            nan_point=find(isnan(d_1));
                            nonnan_point=year;
                            nonnan_point(nan_point)=[];
                            [d_trend,d_detrend]=quadratic_trend(d_1(nonnan_point));
                            tauu_monthly_detrend(nonnan_point,latt,lonn)=d_detrend;
                        else
                            [d_trend,d_detrend]=quadratic_trend(d_1);
                            tauu_monthly_detrend(:,latt,lonn)=d_detrend;
                        end
                    else
                        [d_trend,d_detrend]=quadratic_trend(d);
                        tauu_monthly_detrend(:,latt,lonn)=d_detrend;
                    end
                end
            end
            tauu_monthly=tauu_monthly_detrend;
            clear tauu_monthly_detrend
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tauu_monthly(:,latt,lonn));
                    if numel(find(isnan(d)))>600
                    else
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        d_06=d(6:12:end);
                        d_06=d_06(1:end-1); %%%%%%1900-2098
                        d_07=d(7:12:end);
                        d_07=d_07(1:end-1); %%%%%%1900-2098
                        d_08=d(8:12:end);
                        d_08=d_08(1:end-1); %%%%%%1900-2098
                        d_09=d(9:12:end);
                        d_09=d_09(1:end-1); %%%%%%1900-2098
                        d_10=d(10:12:end);
                        d_10=d_10(1:end-1); %%%%%%1900-2098
                        d_11=d(11:12:end);
                        d_11=d_11(1:end-1); %%%%%%1900-2098
                        d_12=d(12:12:end);
                        d_12=d_12(1:end-1); %%%%%%1900-2098
                        d_01=d(1:12:end);
                        d_01=d_01(2:end); %%%%%%1901-2099
                        d_02=d(2:12:end);
                        d_02=d_02(2:end); %%%%%%1901-2099
                        d_03=d(3:12:end);
                        d_03=d_03(2:end); %%%%%%1901-2099
                        d_04=d(4:12:end);
                        d_04=d_04(2:end); %%%%%%1901-2099
                        d_05=d(5:12:end);
                        d_05=d_05(2:end); %%%%%%1901-2099
                        d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%% JJASONDJFMAM->ONDJF PD
                        d_JJASONDJFMAM_PD=d_JJASONDJFMAM(1:100);
                        tos_nino34_ONDJF_PD=tos_nino34_ONDJF_detrend(num,1:100);
                        tos_nino34_ONDJF_PD=tos_nino34_ONDJF_PD./std_withoutnan(tos_nino34_ONDJF_PD(:));
                        [r rms slope interc n] = ew_regress(tos_nino34_ONDJF_PD,d_JJASONDJFMAM_PD);
                        slope_tauu_J0M1_and_nino34_O0F1_PD(num,latt,lonn)=slope;
                        cor_tauu_J0M1_and_nino34_O0F1_PD(num,latt,lonn)=r;
                        x=d_JJASONDJFMAM_PD;
                        y=tos_nino34_ONDJF_PD;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tauu_J0M1_and_nino34_O0F1_PD(num,latt,lonn)=p;
                        %%%%%%%%%%%%%%%%%%% JJASONDJFMAM->ONDJF, CC
                        d_JJASONDJFMAM_CC=d_JJASONDJFMAM(101:199);
                        tos_nino34_ONDJF_CC=tos_nino34_ONDJF_detrend(num,101:199);
                        tos_nino34_ONDJF_CC=tos_nino34_ONDJF_CC./std_withoutnan(tos_nino34_ONDJF_CC(:));
                        [r rms slope interc n] = ew_regress(tos_nino34_ONDJF_CC,d_JJASONDJFMAM_CC);
                        slope_tauu_J0M1_and_nino34_O0F1_CC(num,latt,lonn)=slope;
                        cor_tauu_J0M1_and_nino34_O0F1_CC(num,latt,lonn)=r;
                        x=d_JJASONDJFMAM_CC;
                        y=tos_nino34_ONDJF_CC;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tauu_J0M1_and_nino34_O0F1_CC(num,latt,lonn)=p;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%% JJASONDJFMAM->DJF PD
                        d_JJASONDJFMAM_PD=d_JJASONDJFMAM(1:100);
                        tos_nino34_DJF_PD=tos_nino34_DJF_detrend(num,1:100);
                        tos_nino34_DJF_PD=tos_nino34_DJF_PD./std_withoutnan(tos_nino34_DJF_PD(:));
                        [r rms slope interc n] = ew_regress(tos_nino34_DJF_PD,d_JJASONDJFMAM_PD);
                        slope_tauu_J0M1_and_nino34_DJF_PD(num,latt,lonn)=slope;
                        cor_tauu_J0M1_and_nino34_DJF_PD(num,latt,lonn)=r;
                        x=d_JJASONDJFMAM_PD;
                        y=tos_nino34_DJF_PD;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tauu_J0M1_and_nino34_DJF_PD(num,latt,lonn)=p;
                        %%%%%%%%%%%%%%%%%%% JJASONDJFMAM->DJF CC
                        d_JJASONDJFMAM_CC=d_JJASONDJFMAM(101:199);
                        tos_nino34_DJF_CC=tos_nino34_DJF_detrend(num,101:199);
                        tos_nino34_DJF_CC=tos_nino34_DJF_CC./std_withoutnan(tos_nino34_DJF_CC(:));
                        [r rms slope interc n] = ew_regress(tos_nino34_DJF_CC,d_JJASONDJFMAM_CC);
                        slope_tauu_J0M1_and_nino34_DJF_CC(num,latt,lonn)=slope;
                        cor_tauu_J0M1_and_nino34_DJF_CC(num,latt,lonn)=r;
                        x=d_JJASONDJFMAM_CC;
                        y=tos_nino34_DJF_CC;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tauu_J0M1_and_nino34_DJF_CC(num,latt,lonn)=p;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear tauu_monthly; 
        end
    end
    save Data_for_Draw_relationship_Tauu_with_Nino34_Seasons_all_CMIP6_New7_corrected_New_New_New.mat models_sst lon_base lat_base ...
         tauu_monthly_mean_s1 tauu_monthly_mean_s2 tauu_djf_mean_s1 tauu_djf_mean_s2 ...
         slope_tauu_J0M1_and_nino34_O0F1_PD cor_tauu_J0M1_and_nino34_O0F1_PD p_tauu_J0M1_and_nino34_O0F1_PD ...   
         slope_tauu_J0M1_and_nino34_DJF_PD  cor_tauu_J0M1_and_nino34_DJF_PD  p_tauu_J0M1_and_nino34_DJF_PD ...
         slope_tauu_J0M1_and_nino34_O0F1_CC cor_tauu_J0M1_and_nino34_O0F1_CC p_tauu_J0M1_and_nino34_O0F1_CC ...
         slope_tauu_J0M1_and_nino34_DJF_CC  cor_tauu_J0M1_and_nino34_DJF_CC  p_tauu_J0M1_and_nino34_DJF_CC 
         
else
    load Data_for_Draw_relationship_Tauu_with_Nino34_Seasons_all_CMIP6_New7_corrected_New_New_New.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load sst_mask.mat
for num=1:numel(models_sst)
    d=squeeze(tauu_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_monthly_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_monthly_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_monthly_mean_s2(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_djf_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_djf_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_djf_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_djf_mean_s2(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(slope_tauu_J0M1_and_nino34_O0F1_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    slope_tauu_J0M1_and_nino34_O0F1_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(cor_tauu_J0M1_and_nino34_O0F1_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    cor_tauu_J0M1_and_nino34_O0F1_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(slope_tauu_J0M1_and_nino34_DJF_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    slope_tauu_J0M1_and_nino34_DJF_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(cor_tauu_J0M1_and_nino34_DJF_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    cor_tauu_J0M1_and_nino34_DJF_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(slope_tauu_J0M1_and_nino34_O0F1_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    slope_tauu_J0M1_and_nino34_O0F1_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(cor_tauu_J0M1_and_nino34_O0F1_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    cor_tauu_J0M1_and_nino34_O0F1_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(slope_tauu_J0M1_and_nino34_DJF_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    slope_tauu_J0M1_and_nino34_DJF_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(cor_tauu_J0M1_and_nino34_DJF_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    cor_tauu_J0M1_and_nino34_DJF_CC(num,:,:)=d;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.2 Calculate sum during EN and LN
if exist('Data_for_Draw_rectification_in_detrendedTauu_v1.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_base=['/datastore/wan143/TAUU/cmip6/monthly_historical_ssp585/1900-2099/tauu_Amon_' cell2mat(models_sst(1)) '_monthly_190001_209912.mat'];
    load(file_base)
    lon_base=lon;
    lat_base=lat;
    clear tauu_monthly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_J0M1_EN0p5_PD =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p5_PD =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_EN0p75_PD=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p75_PD=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_J0M1_EN0p5_CC =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p5_CC =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_EN0p75_CC=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p75_CC=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tauu_J0M1_EN0p5 =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p5 =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_EN0p75=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tauu_J0M1_LN0p75=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/TAUU/cmip6/monthly_historical_ssp585/1900-2099/tauu_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            load(file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(tauu_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Calculate monthly anomalies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tauu_ano =nan(size(tauu_monthly));
            tauu_mean=nan(12,numel(lat),numel(lon));
            for mm=1:12
                tauu_mean(mm,:,:)=squeeze(mean(tauu_monthly(mm:12:1200,:,:),1));
            end
            for tt=1:size(tauu_monthly,1)
                mm=mod(tt,12);
                if mm~=0
                    tauu_ano(tt,:,:)=squeeze(tauu_monthly(tt,:,:))-squeeze(tauu_mean(mm,:,:));
                else
                    tauu_ano(tt,:,:)=squeeze(tauu_monthly(tt,:,:))-squeeze(tauu_mean(12,:,:));
                end
            end
            tauu_monthly=tauu_ano;
            clear tauu_ano
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%Quadratic detrend
            tauu_monthly_detrend=nan(size(tauu_monthly));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tauu_monthly(:,latt,lonn));
                    if numel(find(isnan(d)))>=100
                    elseif numel(find(isnan(d)))>0 & numel(find(isnan(d)))<100
                        d_1=interp1_nan(d);
                        if numel(find(isnan(d_1)))>0
                            year=1:2400;
                            nan_point=find(isnan(d_1));
                            nonnan_point=year;
                            nonnan_point(nan_point)=[];
                            [d_trend,d_detrend]=quadratic_trend(d_1(nonnan_point));
                            tauu_monthly_detrend(nonnan_point,latt,lonn)=d_detrend;
                        else
                            [d_trend,d_detrend]=quadratic_trend(d_1);
                            tauu_monthly_detrend(:,latt,lonn)=d_detrend;
                        end
                    else
                        [d_trend,d_detrend]=quadratic_trend(d);
                        tauu_monthly_detrend(:,latt,lonn)=d_detrend;
                    end
                end
            end
            tauu_monthly=tauu_monthly_detrend;
            clear tauu_monthly_detrend
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tauu_monthly(:,latt,lonn));
                    if numel(find(isnan(d)))>600
                    else
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        d_06=d(6:12:end);
                        d_06=d_06(1:end-1); %%%%%%1900-2098
                        d_07=d(7:12:end);
                        d_07=d_07(1:end-1); %%%%%%1900-2098
                        d_08=d(8:12:end);
                        d_08=d_08(1:end-1); %%%%%%1900-2098
                        d_09=d(9:12:end);
                        d_09=d_09(1:end-1); %%%%%%1900-2098
                        d_10=d(10:12:end);
                        d_10=d_10(1:end-1); %%%%%%1900-2098
                        d_11=d(11:12:end);
                        d_11=d_11(1:end-1); %%%%%%1900-2098
                        d_12=d(12:12:end);
                        d_12=d_12(1:end-1); %%%%%%1900-2098
                        d_01=d(1:12:end);
                        d_01=d_01(2:end); %%%%%%1901-2099
                        d_02=d(2:12:end);
                        d_02=d_02(2:end); %%%%%%1901-2099
                        d_03=d(3:12:end);
                        d_03=d_03(2:end); %%%%%%1901-2099
                        d_04=d(4:12:end);
                        d_04=d_04(2:end); %%%%%%1901-2099
                        d_05=d(5:12:end);
                        d_05=d_05(2:end); %%%%%%1901-2099
                        d_JJASONDJFMAM=(d_06+d_07+d_08+d_09+d_10+d_11+d_12+d_01+d_02+d_03+d_04+d_05)./12;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tauu_J0M1_EN0p5(num,latt,lonn) =sum(d_JJASONDJFMAM(find(tos_nino34_DJF_detrend(num,:)>=0.5)));
                        tauu_J0M1_EN0p75(num,latt,lonn)=sum(d_JJASONDJFMAM(find(tos_nino34_DJF_detrend(num,:)>=0.75)));
                        tauu_J0M1_LN0p5(num,latt,lonn) =sum(d_JJASONDJFMAM(find(tos_nino34_DJF_detrend(num,:)<=-0.5)));
                        tauu_J0M1_LN0p75(num,latt,lonn)=sum(d_JJASONDJFMAM(find(tos_nino34_DJF_detrend(num,:)<=-0.75)));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        d_JJASONDJFMAM_s1           =d_JJASONDJFMAM(1:100);
                        tos_nino34_DJF_detrend_s1   =tos_nino34_DJF_detrend(num,1:100);
                        tauu_J0M1_EN0p5_PD(num,latt,lonn) =sum(d_JJASONDJFMAM_s1(find(tos_nino34_DJF_detrend_s1>=0.5)));
                        tauu_J0M1_EN0p75_PD(num,latt,lonn)=sum(d_JJASONDJFMAM_s1(find(tos_nino34_DJF_detrend_s1>=0.75)));
                        tauu_J0M1_LN0p5_PD(num,latt,lonn) =sum(d_JJASONDJFMAM_s1(find(tos_nino34_DJF_detrend_s1<=-0.5)));
                        tauu_J0M1_LN0p75_PD(num,latt,lonn)=sum(d_JJASONDJFMAM_s1(find(tos_nino34_DJF_detrend_s1<=-0.75)));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        d_JJASONDJFMAM_s2           =d_JJASONDJFMAM(101:end);
                        tos_nino34_DJF_detrend_s2   =tos_nino34_DJF_detrend(num,101:end);
                        tauu_J0M1_EN0p5_CC(num,latt,lonn) =sum(d_JJASONDJFMAM_s2(find(tos_nino34_DJF_detrend_s2>=0.5)));
                        tauu_J0M1_EN0p75_CC(num,latt,lonn)=sum(d_JJASONDJFMAM_s2(find(tos_nino34_DJF_detrend_s2>=0.75)));
                        tauu_J0M1_LN0p5_CC(num,latt,lonn) =sum(d_JJASONDJFMAM_s2(find(tos_nino34_DJF_detrend_s2<=-0.5)));
                        tauu_J0M1_LN0p75_CC(num,latt,lonn)=sum(d_JJASONDJFMAM_s2(find(tos_nino34_DJF_detrend_s2<=-0.75)));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear tauu_monthly;
        end
    end
    save Data_for_Draw_rectification_in_detrendedTauu_v1.mat models_sst lon_base lat_base ... 
         tauu_J0M1_EN0p5 tauu_J0M1_EN0p75 tauu_J0M1_LN0p5 tauu_J0M1_LN0p75 ...
         tauu_J0M1_EN0p5_PD tauu_J0M1_EN0p75_PD tauu_J0M1_LN0p5_PD tauu_J0M1_LN0p75_PD ...
         tauu_J0M1_EN0p5_CC tauu_J0M1_EN0p75_CC tauu_J0M1_LN0p5_CC tauu_J0M1_LN0p75_CC
     
else
    load Data_for_Draw_rectification_in_detrendedTauu_v1.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load sst_mask.mat
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p5(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p5(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p75(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p75(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p5(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p5(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p75(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p75(num,:,:)=d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p5_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p5_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p75_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p75_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p5_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p5_PD(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p75_PD(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p75_PD(num,:,:)=d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p5_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p5_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_EN0p75_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_EN0p75_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p5_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p5_CC(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(tauu_J0M1_LN0p75_CC(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    tauu_J0M1_LN0p75_CC(num,:,:)=d;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Calculate mean state of SST (zonal mean) and regression pattern of PD HUSS onto Nino3.4 (interannual, DJF)
if exist('Data_for_Draw_relationship_SST_with_Nino34_all_CMIP6_New7_corrected.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_base=['/datastore/wan143/tos/CMIP6/monthly_historical_ssp585/1900-2099/tos_Omon_' cell2mat(models_sst(1)) '_monthly_190001_209912.mat'];
    load(file_base)
    lon_base=lon;
    lat_base=lat;
    clear tos_monthly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tos_and_nino34_PD_monthly=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tos_and_nino34_PD_monthly  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tos_and_nino34_PD_monthly    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    slope_tos_and_nino34_PD_djf    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tos_and_nino34_PD_djf      =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tos_and_nino34_PD_djf        =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope_tos_and_nino34_CC_monthly=nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tos_and_nino34_CC_monthly  =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tos_and_nino34_CC_monthly    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    slope_tos_and_nino34_CC_djf    =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    cor_tos_and_nino34_CC_djf      =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    p_tos_and_nino34_CC_djf        =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tos_monthly_mean_s1         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    tos_monthly_mean_s2         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/tos/CMIP6/monthly_historical_ssp585/1900-2099/tos_Omon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            load(file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(tos_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tos_monthly(:,latt,lonn));
                    tos_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    tos_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file=['/datastore/wan143/tos/CMIP6/monthly_historical_ssp585/1900-2099/tos_Omon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file)==0
%             error
        else
            load(file)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(tos_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Calculate monthly anomalies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tos_monthly_ano =nan(size(tos_monthly));
            tos_monthly_mean=nan(12,numel(lat),numel(lon));
            for mm=1:12
                tos_monthly_mean(mm,:,:)=squeeze(mean(tos_monthly(mm:12:1200,:,:),1));
            end
            for tt=1:size(tos_monthly,1)
                mm=mod(tt,12);
                if mm~=0
                    tos_monthly_ano(tt,:,:)=squeeze(tos_monthly(tt,:,:))-squeeze(tos_monthly_mean(mm,:,:));
                else
                    tos_monthly_ano(tt,:,:)=squeeze(tos_monthly(tt,:,:))-squeeze(tos_monthly_mean(12,:,:));
                end
            end
            tos_monthly=tos_monthly_ano;
            clear tos_monthly_ano
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%Quadratic detrend
            tos_monthly_detrend=nan(size(tos_monthly));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tos_monthly(:,latt,lonn));
                    if numel(find(isnan(d)))>=100
                    elseif numel(find(isnan(d)))>0 & numel(find(isnan(d)))<100
                        d_1=interp1_nan(d);
                        if numel(find(isnan(d_1)))>0
                            year=1:2400;
                            nan_point=find(isnan(d_1));
                            nonnan_point=year;
                            nonnan_point(nan_point)=[];
                            [d_trend,d_detrend]=quadratic_trend(d_1(nonnan_point));
                            tos_monthly_detrend(nonnan_point,latt,lonn)=d_detrend;
                        else
                            [d_trend,d_detrend]=quadratic_trend(d_1);
                            tos_monthly_detrend(:,latt,lonn)=d_detrend;
                        end
                    else
                        [d_trend,d_detrend]=quadratic_trend(d);
                        tos_monthly_detrend(:,latt,lonn)=d_detrend;
                    end
                end
            end
            tos_monthly=tos_monthly_detrend;
            clear tos_monthly_detrend
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% Monthly
            tos_nino34_mul_detrend_d=tos_nino34_mul_detrend(num,:);
            tos_nino34_monthly=tos_nino34_mul_detrend_d(1:1200);
            tos_nino34_monthly=tos_nino34_monthly./std_withoutnan(tos_nino34_monthly(:));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tos_monthly(1:1200,latt,lonn));
                    if numel(find(isnan(d)))>600
                    else
                        [r rms slope interc n] = ew_regress(tos_nino34_monthly,d);
                        slope_tos_and_nino34_PD_monthly(num,latt,lonn)=slope;
                        cor_tos_and_nino34_PD_monthly(num,latt,lonn)=r;
                        x=d;
                        y=tos_nino34_monthly;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tos_and_nino34_PD_monthly(num,latt,lonn)  =p;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% Monthly
            tos_nino34_monthly=tos_nino34_mul_detrend_d(1201:end);
            tos_nino34_monthly=tos_nino34_monthly./std_withoutnan(tos_nino34_monthly(:));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(tos_monthly(1201:end,latt,lonn));
                    if numel(find(isnan(d)))>600
                    else
                        [r rms slope interc n] = ew_regress(tos_nino34_monthly,d);
                        slope_tos_and_nino34_CC_monthly(num,latt,lonn)=slope;
                        cor_tos_and_nino34_CC_monthly(num,latt,lonn)=r;
                        x=d;
                        y=tos_nino34_monthly;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tos_and_nino34_CC_monthly(num,latt,lonn)  =p;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% DJF
            tos_nino34_djf_detrend_d=tos_nino34_djf_detrend(num,:);
            tos_nino34_djf=tos_nino34_djf_detrend_d(1:100);
            tos_nino34_djf=tos_nino34_djf./std_withoutnan(tos_nino34_djf(:));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d_d=squeeze(tos_monthly(12:12:end-1,latt,lonn));
                    d_j=squeeze(tos_monthly(13:12:end-1,latt,lonn));
                    d_f=squeeze(tos_monthly(14:12:end-1,latt,lonn));
                    d_djf=(d_d+d_j+d_f)./3;
                    d_djf=d_djf(1:100);
                    if numel(find(isnan(d_djf)))>50
                    else
                        [r rms slope interc n] = ew_regress(tos_nino34_djf,d_djf);
                        slope_tos_and_nino34_PD_djf(num,latt,lonn)=slope;
                        cor_tos_and_nino34_PD_djf(num,latt,lonn)=r;
                        x=d_djf;
                        y=tos_nino34_djf;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tos_and_nino34_PD_djf(num,latt,lonn)  =p;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% DJF
            tos_nino34_djf=tos_nino34_djf_detrend_d(101:end);
            tos_nino34_djf=tos_nino34_djf./std_withoutnan(tos_nino34_djf(:));
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d_d=squeeze(tos_monthly(12:12:end-1,latt,lonn));
                    d_j=squeeze(tos_monthly(13:12:end-1,latt,lonn));
                    d_f=squeeze(tos_monthly(14:12:end-1,latt,lonn));
                    d_djf=(d_d+d_j+d_f)./3;
                    d_djf=d_djf(101:end);
                    if numel(find(isnan(d_djf)))>50
                    else
                        [r rms slope interc n] = ew_regress(tos_nino34_djf,d_djf);
                        slope_tos_and_nino34_CC_djf(num,latt,lonn)=slope;
                        cor_tos_and_nino34_CC_djf(num,latt,lonn)=r;
                        x=d_djf;
                        y=tos_nino34_djf;
                        x(find(isnan(y)))=[];
                        y(find(isnan(y)))=[];
                        y(find(isnan(x)))=[];
                        x(find(isnan(x)))=[];
                        [r p]=my_corr(x(:),y(:));
                        p_tos_and_nino34_CC_djf(num,latt,lonn)  =p;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear tos_monthly; clear tos_monthly_zonal_new
        end
    end
    save Data_for_Draw_relationship_SST_with_Nino34_all_CMIP6_New7_corrected.mat models_sst lon_base lat_base ...
         tos_monthly_mean_s1 tos_monthly_mean_s2 ...
         slope_tos_and_nino34_PD_monthly cor_tos_and_nino34_PD_monthly p_tos_and_nino34_PD_monthly ...
         slope_tos_and_nino34_PD_djf     cor_tos_and_nino34_PD_djf     p_tos_and_nino34_PD_djf ...
         slope_tos_and_nino34_CC_monthly cor_tos_and_nino34_CC_monthly p_tos_and_nino34_CC_monthly ...
         slope_tos_and_nino34_CC_djf     cor_tos_and_nino34_CC_djf     p_tos_and_nino34_CC_djf 
     
else
    load Data_for_Draw_relationship_SST_with_Nino34_all_CMIP6_New7_corrected.mat
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.1 Load in ECS and get rid of models which doesn't have ECS
% load ECS_from_GeraldMeels_ScienceAdvances.mat
% cut_ECS=[];
% i=0;
% for num=1:numel(models_sst)
%     if isnan(ECS(num))==1
%         i=i+1;
%         cut_ECS(i)=num;
%     end
% end
% models_ST(cut_ECS)=[];
% models_sst(cut_ECS)=[];
% ts_GMT_monthly(cut_ECS,:)=[];
% ts_GMT_diff(cut_ECS)     =[];
% ECS(cut_ECS)=[];
% tos_nino34_djf_detrend_std_s1(cut_ECS)=[];
% tos_nino34_djf_detrend_std_s2(cut_ECS)=[];
% thetaoGlobal_monthly_mean_s1(cut_ECS,:,:)=[];
% thetaoGlobal_monthly_mean_s2(cut_ECS,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.1 Load in ECS and get rid of models which doesn't have ECS
cut_tauu_model=[];
i=0;
for num=1:numel(models_sst)
    d=squeeze(tauu_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
        i=i+1;
        cut_tauu_model(i)=num;
    end
end

models_sst(cut_tauu_model) =[];
models_ST(cut_tauu_model)  =[];
ts_GMT_monthly(cut_tauu_model,:)=[];
ts_GMT_diff(cut_tauu_model)     =[];
% ECS(cut_tauu_model)=[];
tos_nino34_djf_detrend_std_s1(cut_tauu_model)=[];
tos_nino34_djf_detrend_std_s2(cut_tauu_model)=[];
tauu_monthly_mean_s1(cut_tauu_model,:,:)=[];
tauu_monthly_mean_s2(cut_tauu_model,:,:)=[];
tauu_djf_mean_s1(cut_tauu_model,:,:)=[];
tauu_djf_mean_s2(cut_tauu_model,:,:)=[];
thetaoGlobal_monthly_mean_s1(cut_tauu_model,:,:)=[];
thetaoGlobal_monthly_mean_s2(cut_tauu_model,:,:)=[];
tos_monthly_mean_s1(cut_tauu_model,:,:)=[];
tos_monthly_mean_s2(cut_tauu_model,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_PD(cut_tauu_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_PD(cut_tauu_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_PD(cut_tauu_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_PD(cut_tauu_model,:,:)=[];     
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_CC(cut_tauu_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_CC(cut_tauu_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_CC(cut_tauu_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_CC(cut_tauu_model,:,:)=[];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_J0M1_EN0p5(cut_tauu_model,:,:)=[];
tauu_J0M1_EN0p75(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p5(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p75(cut_tauu_model,:,:)=[];
tauu_J0M1_EN0p5_PD(cut_tauu_model,:,:)=[];
tauu_J0M1_EN0p75_PD(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p5_PD(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p75_PD(cut_tauu_model,:,:)=[];
tauu_J0M1_EN0p5_CC(cut_tauu_model,:,:)=[];
tauu_J0M1_EN0p75_CC(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p5_CC(cut_tauu_model,:,:)=[];
tauu_J0M1_LN0p75_CC(cut_tauu_model,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.2 Cut off models which is not good {'INM-CM5-0'}    {'MRI-ESM2-0'}
cut_nongood_model=[14 21];

models_sst(cut_nongood_model) =[];
models_ST(cut_nongood_model)  =[];
ts_GMT_monthly(cut_nongood_model,:)=[];
ts_GMT_diff(cut_nongood_model)     =[];
% ECS(cut_nongood_model)=[];
tos_nino34_djf_detrend_std_s1(cut_nongood_model)=[];
tos_nino34_djf_detrend_std_s2(cut_nongood_model)=[];
tauu_monthly_mean_s1(cut_nongood_model,:,:)=[];
tauu_monthly_mean_s2(cut_nongood_model,:,:)=[];
tauu_djf_mean_s1(cut_nongood_model,:,:)=[];
tauu_djf_mean_s2(cut_nongood_model,:,:)=[];
thetaoGlobal_monthly_mean_s1(cut_nongood_model,:,:)=[];
thetaoGlobal_monthly_mean_s2(cut_nongood_model,:,:)=[];
tos_monthly_mean_s1(cut_nongood_model,:,:)=[];
tos_monthly_mean_s2(cut_nongood_model,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_PD(cut_nongood_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_PD(cut_nongood_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_PD(cut_nongood_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_PD(cut_nongood_model,:,:)=[];     
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_CC(cut_nongood_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_CC(cut_nongood_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_CC(cut_nongood_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_CC(cut_nongood_model,:,:)=[];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_J0M1_EN0p5(cut_nongood_model,:,:)=[];
tauu_J0M1_EN0p75(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p5(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p75(cut_nongood_model,:,:)=[];
tauu_J0M1_EN0p5_PD(cut_nongood_model,:,:)=[];
tauu_J0M1_EN0p75_PD(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p5_PD(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p75_PD(cut_nongood_model,:,:)=[];
tauu_J0M1_EN0p5_CC(cut_nongood_model,:,:)=[];
tauu_J0M1_EN0p75_CC(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p5_CC(cut_nongood_model,:,:)=[];
tauu_J0M1_LN0p75_CC(cut_nongood_model,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.3 Load in ECS and get rid of models which doesn't have thetao
cut_thetao_model=[];
i=0;
for num=1:numel(models_sst)
    d=squeeze(thetaoGlobal_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
        i=i+1;
        cut_thetao_model(i)=num;
    end
end

models_sst(cut_thetao_model) =[];
models_ST(cut_thetao_model)  =[];
ts_GMT_monthly(cut_thetao_model,:)=[];
ts_GMT_diff(cut_thetao_model)     =[];
% ECS(cut_thetao)=[];
tos_nino34_djf_detrend_std_s1(cut_thetao_model)=[];
tos_nino34_djf_detrend_std_s2(cut_thetao_model)=[];
tauu_monthly_mean_s1(cut_thetao_model,:,:)=[];
tauu_monthly_mean_s2(cut_thetao_model,:,:)=[];
tauu_djf_mean_s1(cut_thetao_model,:,:)=[];
tauu_djf_mean_s2(cut_thetao_model,:,:)=[];
thetaoGlobal_monthly_mean_s1(cut_thetao_model,:,:)=[];
thetaoGlobal_monthly_mean_s2(cut_thetao_model,:,:)=[];
tos_monthly_mean_s1(cut_thetao_model,:,:)=[];
tos_monthly_mean_s2(cut_thetao_model,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_PD(cut_thetao_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_PD(cut_thetao_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_PD(cut_thetao_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_PD(cut_thetao_model,:,:)=[];     
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_CC(cut_thetao_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_CC(cut_thetao_model,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_CC(cut_thetao_model,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_CC(cut_thetao_model,:,:)=[];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_J0M1_EN0p5(cut_thetao_model,:,:)=[];
tauu_J0M1_EN0p75(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p5(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p75(cut_thetao_model,:,:)=[];
tauu_J0M1_EN0p5_PD(cut_thetao_model,:,:)=[];
tauu_J0M1_EN0p75_PD(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p5_PD(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p75_PD(cut_thetao_model,:,:)=[];
tauu_J0M1_EN0p5_CC(cut_thetao_model,:,:)=[];
tauu_J0M1_EN0p75_CC(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p5_CC(cut_thetao_model,:,:)=[];
tauu_J0M1_LN0p75_CC(cut_thetao_model,:,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Projected changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tos_nino34_djf_detrend_std_diff=tos_nino34_djf_detrend_std_s2-tos_nino34_djf_detrend_std_s1;
for num=1:numel(models_sst)
    tos_nino34_djf_detrend_std_diff(num)=tos_nino34_djf_detrend_std_diff(num)./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_monthly_mean_diff=tauu_monthly_mean_s2-tauu_monthly_mean_s1;
for num=1:numel(models_sst)
    tauu_monthly_mean_diff(num,:,:)=squeeze(tauu_monthly_mean_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_djf_mean_diff=tauu_djf_mean_s2-tauu_djf_mean_s1;
for num=1:numel(models_sst)
    tauu_djf_mean_diff(num,:,:)=squeeze(tauu_djf_mean_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaoGlobal_monthly_mean_diff=thetaoGlobal_monthly_mean_s2-thetaoGlobal_monthly_mean_s1;
for num=1:numel(models_sst)
    thetaoGlobal_monthly_mean_diff(num,:,:)=squeeze(thetaoGlobal_monthly_mean_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tos_monthly_mean_diff=tos_monthly_mean_s2-tos_monthly_mean_s1;
for num=1:numel(models_sst)
    tos_monthly_mean_diff(num,:,:)=squeeze(tos_monthly_mean_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_diff=slope_tauu_J0M1_and_nino34_O0F1_CC-slope_tauu_J0M1_and_nino34_O0F1_PD;
for num=1:numel(models_sst)
    slope_tauu_J0M1_and_nino34_O0F1_diff(num,:,:)=squeeze(slope_tauu_J0M1_and_nino34_O0F1_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cor_tauu_J0M1_and_nino34_O0F1_diff=cor_tauu_J0M1_and_nino34_O0F1_CC-cor_tauu_J0M1_and_nino34_O0F1_PD;
for num=1:numel(models_sst)
    cor_tauu_J0M1_and_nino34_O0F1_diff(num,:,:)=squeeze(cor_tauu_J0M1_and_nino34_O0F1_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_DJF_diff=slope_tauu_J0M1_and_nino34_DJF_CC-slope_tauu_J0M1_and_nino34_DJF_PD;
for num=1:numel(models_sst)
    slope_tauu_J0M1_and_nino34_DJF_diff(num,:,:)=squeeze(slope_tauu_J0M1_and_nino34_DJF_diff(num,:,:))./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cor_tauu_J0M1_and_nino34_DJF_diff=cor_tauu_J0M1_and_nino34_DJF_CC-cor_tauu_J0M1_and_nino34_DJF_PD;
for num=1:numel(models_sst)
    cor_tauu_J0M1_and_nino34_DJF_diff(num,:,:)=squeeze(cor_tauu_J0M1_and_nino34_DJF_diff(num,:,:))./ts_GMT_diff(num);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Calculate SOI
lev_deeper_region=1000;
lev_shallower_region=10;
lat_lower_region=-60;
lat_upper_region=-40;
thetaoGlobal_monthly_mean_diff_region=nan(numel(models_sst),1);
for num=1:numel(models_sst)
    d=squeeze(thetaoGlobal_monthly_mean_diff(num,find(lev_thetao>=lev_shallower_region & lev_thetao<=lev_deeper_region),find(lat_thetao>=lat_lower_region & lat_thetao<=lat_upper_region)));
    thetaoGlobal_monthly_mean_diff_region(num)=mean_withoutnan(d(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.1 Draw regression pattern of tauu onto Nino3.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_tauu_onto_Nino34=nan(numel(lat_base),numel(lon_base));
cor_tauu_onto_Nino34=nan(numel(lat_base),numel(lon_base));
p_tauu_onto_Nino34  =nan(numel(lat_base),numel(lon_base));

for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_monthly_mean_diff(:,latt,lonn));
        if numel(find(isnan(d)))>=numel(d)*0.5
        else
            [r rms slope interc n] = ew_regress(tos_nino34_djf_detrend_std_diff,d);
            reg_tauu_onto_Nino34(latt,lonn)=slope;
            x=tos_nino34_djf_detrend_std_diff;
            x(find(isnan(d)))=[];
            d(find(isnan(d)))=[];
            [r p]=my_corr(x(:),d(:));
            cor_tauu_onto_Nino34(latt,lonn)=r;
            p_tauu_onto_Nino34(latt,lonn)  =p;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'position',[77          25        1521         971])
load mmjet.mat

lat_lower=-80;
lat_upper=-20;
lon_left=0;
lon_right=360;

lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));
lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));

subplot(2,2,1)
reg_tauu_onto_Nino34_d=reg_tauu_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
[c h]=contourf(lon_1,lat_1,reg_tauu_onto_Nino34_d,[-1:0.001:1.0]);
set(h,'color','none')
caxis([-0.021 0.021])
colorbar('Ticks',[-0.018 -0.012 -0.006 0 0.006 0.012 0.018])
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
text(lon_left,lat_upper+5,'{\bfA} Changes in westerlies due to changes in ENSO amplitude','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
set(gca,'position',[0.18 0.7192 0.3222    0.2357])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_tauu_onto_Nino34_d=p_tauu_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
for latt=1:2:numel(lat_1)
    for lonn=1:2:numel(lon_1)
        if p_tauu_onto_Nino34_d(latt,lonn)<0.1
            plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.2 Compare top5 and bottom 5, rectification.
lat_lower=-80;
lat_upper=-20;
lon_left=0;
lon_right=360;

lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));
lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% 0.75
tauu_J0M1_total0p75_PD=tauu_J0M1_EN0p75_PD+tauu_J0M1_LN0p75_PD;
tauu_J0M1_total0p75_CC=tauu_J0M1_EN0p75_CC+tauu_J0M1_LN0p75_CC;
tauu_J0M1_total0p75_diff=tauu_J0M1_total0p75_CC-tauu_J0M1_total0p75_PD;
for num=1:numel(models_sst)
    tauu_J0M1_total0p75_diff(num,:,:)=tauu_J0M1_total0p75_diff(num,:,:)./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%% Bottom 5
[B I]=sort(tos_nino34_djf_detrend_std_diff);
tauu_J0M1_total0p75_diff_bottom5=tauu_J0M1_total0p75_diff(I(1:5),:,:);
tauu_J0M1_total0p75_diff_bottom5_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p75_diff_bottom5(:,latt,lonn));
        tauu_J0M1_total0p75_diff_bottom5_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p75_diff_bottom5_mul_d=tauu_J0M1_total0p75_diff_bottom5_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p75_diff_bottom5_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p75_diff_bottom5_mul_d(latt,:);
    tauu_J0M1_total0p75_diff_bottom5_mul_d_zonal(latt)=mean_withoutnan(d);
end
%%%%%%%%%%%%%%%%%%%%%% Top 5
tauu_J0M1_total0p75_diff_top5=tauu_J0M1_total0p75_diff(I(end-4:end),:,:);
tauu_J0M1_total0p75_diff_top5_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p75_diff_top5(:,latt,lonn));
        tauu_J0M1_total0p75_diff_top5_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p75_diff_top5_mul_d=tauu_J0M1_total0p75_diff_top5_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p75_diff_top5_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p75_diff_top5_mul_d(latt,:);
    tauu_J0M1_total0p75_diff_top5_mul_d_zonal(latt)=mean_withoutnan(d);
end
%%%%%%%%%%%%%%%%%%%%%% Difference
tauu_J0M1_total0p75_diff_TminusB=tauu_J0M1_total0p75_diff_top5-tauu_J0M1_total0p75_diff_bottom5;
tauu_J0M1_total0p75_diff_TminusB_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p75_diff_TminusB(:,latt,lonn));
        tauu_J0M1_total0p75_diff_TminusB_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p75_diff_TminusB_mul_d=tauu_J0M1_total0p75_diff_TminusB_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p75_diff_TminusB_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p75_diff_TminusB_mul_d(latt,:);
    tauu_J0M1_total0p75_diff_TminusB_mul_d_zonal(latt)=mean_withoutnan(d);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% 0.5
tauu_J0M1_total0p5_PD=tauu_J0M1_EN0p5_PD+tauu_J0M1_LN0p5_PD;
tauu_J0M1_total0p5_CC=tauu_J0M1_EN0p5_CC+tauu_J0M1_LN0p5_CC;
tauu_J0M1_total0p5_diff=tauu_J0M1_total0p5_CC-tauu_J0M1_total0p5_PD;
for num=1:numel(models_sst)
    tauu_J0M1_total0p5_diff(num,:,:)=tauu_J0M1_total0p5_diff(num,:,:)./ts_GMT_diff(num);
end
%%%%%%%%%%%%%%%%%%%%%% Bottom 5
[B I]=sort(tos_nino34_djf_detrend_std_diff);
tauu_J0M1_total0p5_diff_bottom5=tauu_J0M1_total0p5_diff(I(1:5),:,:);
tauu_J0M1_total0p5_diff_bottom5_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p5_diff_bottom5(:,latt,lonn));
        tauu_J0M1_total0p5_diff_bottom5_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p5_diff_bottom5_mul_d=tauu_J0M1_total0p5_diff_bottom5_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p5_diff_bottom5_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p5_diff_bottom5_mul_d(latt,:);
    tauu_J0M1_total0p5_diff_bottom5_mul_d_zonal(latt)=mean_withoutnan(d);
end
%%%%%%%%%%%%%%%%%%%%%% Top 5
tauu_J0M1_total0p5_diff_top5=tauu_J0M1_total0p5_diff(I(end-4:end),:,:);
tauu_J0M1_total0p5_diff_top5_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p5_diff_top5(:,latt,lonn));
        tauu_J0M1_total0p5_diff_top5_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p5_diff_top5_mul_d=tauu_J0M1_total0p5_diff_top5_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p5_diff_top5_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p5_diff_top5_mul_d(latt,:);
    tauu_J0M1_total0p5_diff_top5_mul_d_zonal(latt)=mean_withoutnan(d);
end
%%%%%%%%%%%%%%%%%%%%%% Difference
tauu_J0M1_total0p5_diff_TminusB=tauu_J0M1_total0p5_diff_top5-tauu_J0M1_total0p5_diff_bottom5;
tauu_J0M1_total0p5_diff_TminusB_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_J0M1_total0p5_diff_TminusB(:,latt,lonn));
        tauu_J0M1_total0p5_diff_TminusB_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_J0M1_total0p5_diff_TminusB_mul_d=tauu_J0M1_total0p5_diff_TminusB_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
tauu_J0M1_total0p5_diff_TminusB_mul_d_zonal=nan(numel(lat_1),1);
for latt=1:numel(lat_1)
    d=tauu_J0M1_total0p5_diff_TminusB_mul_d(latt,:);
    tauu_J0M1_total0p5_diff_TminusB_mul_d_zonal(latt)=mean_withoutnan(d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
plot(tauu_J0M1_total0p75_diff_top5_mul_d_zonal,lat_1,'color','r','linestyle','-','linewidth',2)
hold on
plot(tauu_J0M1_total0p75_diff_bottom5_mul_d_zonal,lat_1,'color','b','linestyle','-','linewidth',2)
plot(tauu_J0M1_total0p75_diff_TminusB_mul_d_zonal,lat_1,'color','k','linestyle','-','linewidth',3)
% plot(tauu_J0M1_total0p5_diff_TminusB_mul_d_zonal,lat_1,'color','k','linestyle','--','linewidth',3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'tickdir','out')
set(gca,'fontsize',15)
box on
text(-0.06,-15,'{\bfB} Changes in cumulative ENSO anomalies','fontsize',15)
xlim([-0.06 0.06])
set(gca,'xtick',[-0.06:0.02:0.06],'xticklabel',[-0.06:0.02:0.06])
ylim([lat_lower lat_upper])
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
set(gca,'xminortick','on','yminortick','on')
hold on
plot(zeros(numel(lat_1),1),lat_1,'k--','linewidth',3)
% set(gca,'position',[0.5400    0.5838    0.2879    0.3412])
set(gca,'position',[0.18 0.0917 0.3222    0.4367])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fill_region=-0.1:0.01:0.1;
fill([fill_region,fliplr(fill_region)],[ones(1,numel(fill_region)).*-60,fliplr(ones(1,numel(fill_region)).*-40)],[0.85 0.85 0.85],'LineStyle','none','FaceAlpha',0.5)
legend('Top 5','Bottom 5','Difference')

































