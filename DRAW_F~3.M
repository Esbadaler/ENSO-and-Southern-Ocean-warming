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
%% 3. Calculate mean state of Hur (zonal mean) and regression pattern of PD Hur onto Nino3.4 (interannual, DJF)
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
%% 4.1 Load in ECS and get rid of models which doesn't have ECS
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.2 Cut off models which is not good {'INM-CM5-0'}    {'MRI-ESM2-0'}
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.3 Load in ECS and get rid of models which doesn't have thetao
cut_thetao=[];
i=0;
for num=1:numel(models_sst)
    d=squeeze(thetaoGlobal_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
        i=i+1;
        cut_thetao(i)=num;
    end
end

models_sst(cut_thetao) =[];
models_ST(cut_thetao)  =[];
ts_GMT_monthly(cut_thetao,:)=[];
ts_GMT_diff(cut_thetao)     =[];
% ECS(cut_thetao)=[];
tos_nino34_djf_detrend_std_s1(cut_thetao)=[];
tos_nino34_djf_detrend_std_s2(cut_thetao)=[];
tauu_monthly_mean_s1(cut_thetao,:,:)=[];
tauu_monthly_mean_s2(cut_thetao,:,:)=[];
tauu_djf_mean_s1(cut_thetao,:,:)=[];
tauu_djf_mean_s2(cut_thetao,:,:)=[];
thetaoGlobal_monthly_mean_s1(cut_thetao,:,:)=[];
thetaoGlobal_monthly_mean_s2(cut_thetao,:,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_PD(cut_thetao,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_PD(cut_thetao,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_PD(cut_thetao,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_PD(cut_thetao,:,:)=[];     
%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope_tauu_J0M1_and_nino34_O0F1_CC(cut_thetao,:,:)=[];   
cor_tauu_J0M1_and_nino34_O0F1_CC(cut_thetao,:,:)=[];     
slope_tauu_J0M1_and_nino34_DJF_CC(cut_thetao,:,:)=[];   
cor_tauu_J0M1_and_nino34_DJF_CC(cut_thetao,:,:)=[];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Projected changes
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
%% 6. Calculate SOI
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
%% 7. EOF in {\Delta} ENSO-tauu
load mmjet.mat

kmod=3;
lat_lower=-80;
lat_upper=-20;
lon_left=0;
lon_right=360;

lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));
lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));

% cut_model=[10 11];
cut_model=[];
slope_tauu_J0M1_and_nino34_DJF_diff(cut_model,:,:)=[];
tos_nino34_djf_detrend_std_diff(cut_model)=[];
models_new=models_sst;
models_new(cut_model)=[];

slope_tauu_J0M1_and_nino34_DJF_diff_d=slope_tauu_J0M1_and_nino34_DJF_diff(:,find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
[V_eof,PC,per,lamd]=EOF_after_detrended(slope_tauu_J0M1_and_nino34_DJF_diff_d,kmod);

pc_1=PC(:,1);
pc_2=PC(:,2);
pattern_1=squeeze(V_eof(:,:,1));
pattern_2=squeeze(V_eof(:,:,2));
pattern_1=pattern_1.*std(pc_1);
pc_1     =pc_1./std(pc_1);
pattern_2=pattern_2.*std(pc_2);
pc_2     =pc_2./std(pc_2);

[r rms slope interc n] = ew_regress(pc_1,tos_nino34_djf_detrend_std_diff(:));
r
[r rms slope interc n] = ew_regress(pc_2,tos_nino34_djf_detrend_std_diff(:));
r

figure
set(gcf,'position',[77          25        1521         971])

subplot(2,2,1)
[c h]=contourf(lon_1,lat_1,pattern_1,[-0.1:0.0001:0.1]);
set(h,'color','none')
caxis([-0.014 0.014].*0.1)
colorbar('Ticks',[-0.012 -0.008 -0.004 0 0.004 0.008 0.012].*0.1)
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
%set(gca,'ydir','reverse')
set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
text(lon_left,lat_upper+5,'{\bfA} Pacific-South American Pattern','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
set(gca,'position',[0.18 0.7192 0.3222    0.2357])

subplot(2,2,3)
h=nan(numel(models_new),1);
for num=1:numel(models_new)
    [mar,fco,si,edge]=ew_CMIP_markers_wang(num);
    h(num)=plot(pc_1(num),thetaoGlobal_monthly_mean_diff_region(num),'marker',mar,'markerfacecolor',fco,'markersize',si+6,'markeredgecolor',edge,'linestyle','none');
    hold on
end
xlabel('Inter-model difference in PSA pattern')
ylabel('SO warming ^oC per degree C of GW')
xlim([-2 4])
ylim([0.2 0.65])
set(gca,'xtick',[-2:1:3],'xticklabel',[-2:1:3])
set(gca,'ytick',[0.2:0.1:0.6],'yticklabel',[0.2:0.1:0.6])
text(-2,0.68,'{\bfB} PSA vs SO warming','fontsize',15)
axis normal
set(gca,'tickdir','out')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r rms slope interc n] = ew_regress(pc_1,thetaoGlobal_monthly_mean_diff_region);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(pc_1,pc_1.*slope+interc,'k','linewidth',2)
x=nan(numel(pc_1)+2,1);
y=nan(numel(thetaoGlobal_monthly_mean_diff_region)+2,1);
x(2:end-1)=pc_1;
y(2:end-1)=thetaoGlobal_monthly_mean_diff_region;
x(1)      =-2;
y(1)      =-2*slope+interc;
y(end)    =0.2;
x(end)    =(0.2-interc)/slope;
plot(x,x.*slope+interc,'k','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text(-1.9,0.33,['Corre. Coeff.= ' num2str(r,2)],'fontsize',15)
text(-1.9,0.28,['Slope= ' num2str(slope,2)],'fontsize',15)
[r p]=my_corr(pc_1(:),thetaoGlobal_monthly_mean_diff_region(:));
if p<0.001
    text(-1.9,0.23,'p<0.001','fontsize',15)
else
    text(-1.9,0.23,['p= ' num2str(p,2)],'fontsize',15)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'fontsize',15)
set(gca,'xminortick','on','yminortick','on')
plot(-100:1:100,ones(numel(-100:1:100),1).*0,'--k')
plot(ones(numel(-100:1:100),1).*0,-100:1:100,'--k')
set(gca,'position',[0.18 0.0917 0.3222    0.4367])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=pc_1;
y=thetaoGlobal_monthly_mean_diff_region;

out_vec=fit_ellipse_new(x,y);

z=[out_vec(1) out_vec(2)];
plotellipse(z, out_vec(3), out_vec(4), out_vec(5))
legend(h,models_new)



figure
set(gcf,'position',[77          25        1521         971])
subplot(2,2,1)
[c h]=contourf(lon_1,lat_1,pattern_2,[-0.1:0.0001:0.1]);
set(h,'color','none')
caxis([-0.028 0.028].*0.1)
colorbar('Ticks',[-0.024 -0.016 -0.008 0 0.008 0.016 0.024].*0.1)
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
%set(gca,'ydir','reverse')
set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
text(lon_left,lat_upper+5,'{\bfA} Zonally symmetric pattern','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lon_left_region=0;
% lon_right_region=360;
% lat_lower_region=-70;
% lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)
% pattern_region=pattern_2(find(lat_1>=lat_lower_region & lat_1<=lat_upper_region),find(lon_1>=lon_left_region & lon_1<=lon_right_region));
% pattern_region=mean_withoutnan(pattern_region(:));
set(gca,'position',[0.18 0.7192 0.3222    0.2357])


subplot(2,2,3)
h=nan(numel(models_new),1);
for num=1:numel(models_new)
    [mar,fco,si,edge]=ew_CMIP_markers_wang(num);
    h(num)=plot(pc_2(num),thetaoGlobal_monthly_mean_diff_region(num),'marker',mar,'markerfacecolor',fco,'markersize',si+6,'markeredgecolor',edge,'linestyle','none');
    hold on
end
xlabel('Inter-model difference in zonally symmetric pattern')
ylabel('SO warming ^oC per degree C of GW')
xlim([-4 3])
ylim([0.2 0.65])
set(gca,'xtick',[-4:1:3],'xticklabel',[-4:1:3])
set(gca,'ytick',[0.2:0.1:0.6],'yticklabel',[0.2:0.1:0.6])
text(-4,0.68,'{\bfB} Zonally symmetric pattern vs SO warming','fontsize',15)
axis normal
set(gca,'tickdir','out')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r rms slope interc n] = ew_regress(pc_2,thetaoGlobal_monthly_mean_diff_region);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(pc_1,pc_1.*slope+interc,'k','linewidth',2)
x=nan(numel(pc_2)+2,1);
y=nan(numel(thetaoGlobal_monthly_mean_diff_region)+2,1);
x(2:end-1)=pc_2;
y(2:end-1)=thetaoGlobal_monthly_mean_diff_region;
x(1)      =-4;
y(1)      =-4*slope+interc;
y(end)    =0.2;
x(end)    =(0.2-interc)/slope;
plot(x,x.*slope+interc,'k','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(pc_2,pc_2.*slope+interc,'k','linewidth',2)
text(-3.9,0.33,['Corre. Coeff.= ' num2str(r,2)],'fontsize',15)
text(-3.9,0.28,['Slope= ' num2str(slope,2)],'fontsize',15)
[r p]=my_corr(pc_2(:),thetaoGlobal_monthly_mean_diff_region(:));
if p<0.001
    text(-3.9,0.23,'p<0.001','fontsize',15)
else
    text(-3.9,0.23,['p= ' num2str(p,2)],'fontsize',15)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'fontsize',15)
set(gca,'xminortick','on','yminortick','on')
plot(-100:1:100,ones(numel(-100:1:100),1).*0,'--k')
plot(ones(numel(-100:1:100),1).*0,-100:1:100,'--k')
set(gca,'position',[0.18 0.0917 0.3222    0.4367])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=pc_2;
y=thetaoGlobal_monthly_mean_diff_region;

out_vec=fit_ellipse_new(x,y);

z=[out_vec(1) out_vec(2)];
plotellipse(z, out_vec(3), out_vec(4), out_vec(5))
legend(h,models_sst)



























































