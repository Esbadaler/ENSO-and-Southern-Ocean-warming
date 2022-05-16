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
%% 3. Calculate mean state of SW, LW, LHF, and SHF
if exist('Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_base=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/rsds/monthly_historical_ssp585/1900-2099/rsds_Amon_' cell2mat(models_sst(1)) '_monthly_190001_209912.mat'];
    load(file_base)
    lon_base=lon;
    lat_base=lat;
    clear rsds_monthly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sw_monthly_mean_s1         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    sw_monthly_mean_s2         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file_rsds=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/rsds/monthly_historical_ssp585/1900-2099/rsds_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        file_rsus=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/rsus/monthly_historical_ssp585/1900-2099/rsus_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file_rsds)==0 | exist(file_rsus)==0
%             error
        else
            load(file_rsds)
            load(file_rsus)
            sw_monthly=rsds_monthly-rsus_monthly;
            clear rsds_monthly; clear rsus_monthly
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(sw_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(sw_monthly(:,latt,lonn));
                    sw_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    sw_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear sw_monthly; clear file_rsds; clear file_rsus
        end
    end
    save Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat models_sst lon_base lat_base ...
         sw_monthly_mean_s1 sw_monthly_mean_s2 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lw_monthly_mean_s1         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    lw_monthly_mean_s2         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file_rlds=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/rlds/monthly_historical_ssp585/1900-2099/rlds_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        file_rlus=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/rlus/monthly_historical_ssp585/1900-2099/rlus_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file_rlds)==0 | exist(file_rlus)==0
%             error
        else
            load(file_rlds)
            load(file_rlus)
            lw_monthly=rlds_monthly-rlus_monthly;
            clear rlds_monthly; clear rlus_monthly
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(lw_monthly,1)~=2400
                error
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(lw_monthly(:,latt,lonn));
                    lw_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    lw_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear lw_monthly; clear file_rlds; clear file_rlus
        end
    end
    save Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat models_sst lon_base lat_base ...
         lw_monthly_mean_s1 lw_monthly_mean_s2  -append
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfls_monthly_mean_s1         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    hfls_monthly_mean_s2         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file_hfls=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/hfls/monthly_historical_ssp585/1900-2099/hfls_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file_hfls)==0
%             error
        else
            load(file_hfls)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(hfls_monthly,1)~=2400
                error
            end
            hfls_monthly=hfls_monthly.*-1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(hfls_monthly(:,latt,lonn));
                    hfls_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    hfls_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear hfls_monthly; clear file_hfls
        end
    end
    save Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat models_sst lon_base lat_base ...
         hfls_monthly_mean_s1 hfls_monthly_mean_s2  -append
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfss_monthly_mean_s1         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    hfss_monthly_mean_s2         =nan(numel(models_sst),numel(lat_base),numel(lon_base));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num=1:numel(models_sst)
        num
        file_hfss=['/datastore/wan143/Ocean_data_related_to_HeatBudget/CMIP6/hfss/monthly_historical_ssp585/1900-2099/hfss_Amon_' cell2mat(models_sst(num)) '_monthly_190001_209912.mat'];
        if exist(file_hfss)==0
%             error
        else
            load(file_hfss)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(hfss_monthly,1)~=2400
                error
            end
            hfss_monthly=hfss_monthly.*-1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for latt=1:numel(lat_base)
                for lonn=1:numel(lon_base)
                    d=squeeze(hfss_monthly(:,latt,lonn));
                    hfss_monthly_mean_s1(num,latt,lonn)=mean_withoutnan(d(1:1200));
                    hfss_monthly_mean_s2(num,latt,lonn)=mean_withoutnan(d(1201:end));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear hfss_monthly; clear file_hfss
        end
    end
    save Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat models_sst lon_base lat_base ...
         hfss_monthly_mean_s1 hfss_monthly_mean_s2  -append
     
else
    load Draw_relationship_wap_with_Nino34_all_CMIP6_New7_corrected_Final_with_Qnet.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load sst_mask.mat
for num=1:numel(models_sst)
    d=squeeze(sw_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    sw_monthly_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(sw_monthly_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    sw_monthly_mean_s2(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(lw_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    lw_monthly_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(lw_monthly_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    lw_monthly_mean_s2(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(hfls_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    hfls_monthly_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(hfls_monthly_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    hfls_monthly_mean_s2(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(hfss_monthly_mean_s1(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    hfss_monthly_mean_s1(num,:,:)=d;
end
for num=1:numel(models_sst)
    d=squeeze(hfss_monthly_mean_s2(num,:,:));
    if numel(find(isnan(d)))==numel(d)
    else
        d(find(isnan(sst_mask)==1))=nan;
    end
    hfss_monthly_mean_s2(num,:,:)=d;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Calculate mean state of Hur (zonal mean) and regression pattern of PD Hur onto Nino3.4 (interannual, DJF)
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
sw_monthly_mean_s1(cut_tauu_model,:,:)=[];
sw_monthly_mean_s2(cut_tauu_model,:,:)=[];
lw_monthly_mean_s1(cut_tauu_model,:,:)=[];
lw_monthly_mean_s2(cut_tauu_model,:,:)=[];
hfls_monthly_mean_s1(cut_tauu_model,:,:)=[];
hfls_monthly_mean_s2(cut_tauu_model,:,:)=[];
hfss_monthly_mean_s1(cut_tauu_model,:,:)=[];
hfss_monthly_mean_s2(cut_tauu_model,:,:)=[];

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
sw_monthly_mean_s1(cut_nongood_model,:,:)=[];
sw_monthly_mean_s2(cut_nongood_model,:,:)=[];
lw_monthly_mean_s1(cut_nongood_model,:,:)=[];
lw_monthly_mean_s2(cut_nongood_model,:,:)=[];
hfls_monthly_mean_s1(cut_nongood_model,:,:)=[];
hfls_monthly_mean_s2(cut_nongood_model,:,:)=[];
hfss_monthly_mean_s1(cut_nongood_model,:,:)=[];
hfss_monthly_mean_s2(cut_nongood_model,:,:)=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5.3 Load in ECS and get rid of models which doesn't have thetao
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
sw_monthly_mean_s1(cut_thetao,:,:)=[];
sw_monthly_mean_s2(cut_thetao,:,:)=[];
lw_monthly_mean_s1(cut_thetao,:,:)=[];
lw_monthly_mean_s2(cut_thetao,:,:)=[];
hfls_monthly_mean_s1(cut_thetao,:,:)=[];
hfls_monthly_mean_s2(cut_thetao,:,:)=[];
hfss_monthly_mean_s1(cut_thetao,:,:)=[];
hfss_monthly_mean_s2(cut_thetao,:,:)=[];

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
Qnet_monthly_mean_s1=sw_monthly_mean_s1+lw_monthly_mean_s1+hfls_monthly_mean_s1+hfss_monthly_mean_s1;
Qnet_monthly_mean_s2=sw_monthly_mean_s2+lw_monthly_mean_s2+hfls_monthly_mean_s2+hfss_monthly_mean_s2;

Qnet_monthly_mean_diff=Qnet_monthly_mean_s2-Qnet_monthly_mean_s1;
for num=1:numel(models_sst)
    Qnet_monthly_mean_diff(num,:,:)=squeeze(Qnet_monthly_mean_diff(num,:,:))./ts_GMT_diff(num);
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
%% 8.1 Projected changes in thetao, scaled.
thetaoGlobal_monthly_mean_diff=thetaoGlobal_monthly_mean_s2-thetaoGlobal_monthly_mean_s1;
for num=1:numel(models_sst)
    thetaoGlobal_monthly_mean_diff(num,:,:)=squeeze(thetaoGlobal_monthly_mean_diff(num,:,:))./ts_GMT_diff(num);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lev_deeper=2000;
lev_shallower=10;
lat_lower=-72;
lat_upper=-18;

lev_1=lev_thetao(find(lev_thetao>=lev_shallower & lev_thetao<=lev_deeper));
lat_1=lat_thetao(find(lat_thetao>=lat_lower     & lat_thetao<=lat_upper));

load mmjet.mat
figure
set(gcf,'position',[1          31        1920         973])

subplot(2,2,3)
thetaoGlobal_monthly_mean_diff_ave=nan(numel(lev_thetao),numel(lat_thetao));
for levv=1:numel(lev_thetao)
    for latt=1:numel(lat_thetao)
        d=squeeze(thetaoGlobal_monthly_mean_diff(:,levv,latt));
        thetaoGlobal_monthly_mean_diff_ave(levv,latt)=mean_withoutnan(d);
    end
end
thetaoGlobal_monthly_mean_diff_ave_d=thetaoGlobal_monthly_mean_diff_ave(find(lev_thetao>=lev_shallower & lev_thetao<=lev_deeper),find(lat_thetao>=lat_lower & lat_thetao<=lat_upper));
[c h]=contourf(lat_1,lev_1,thetaoGlobal_monthly_mean_diff_ave_d,[-1:0.001:1.0]);
set(h,'color','none')
caxis([0 0.7])
colorbar('Ticks',[0 0.2 0.4 0.6])
colormap(gca,mm_jet_normal(8:end,:))
set(gca,'tickdir','out')
set(gca,'ydir','reverse')
ylabel('m')
% set(gca,'xtick',[-80:10:80],'xticklabel',{'80^oS','70^oS','60^oS','50^oS','40^oS','30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN','80^oN'})
% set(gca,'ytick',[10 50 100 200 300 400 500 600 700 800 900],'yticklabel',[10 50 100 200 300 400 500 600 700 800 900])
set(gca,'xtick',[-70:5:-20],'xticklabel',{'70^oS','','60^oS','','50^oS','','40^oS','','30^oS','','20^oS'})
set(gca,'ytick',[100:200:2000],'yticklabel',[100:200:2000])
text(-70,-100,'{\bfb} Multi-model ensemble mean temperature changes','fontsize',15)
hold on
xlim([-70 -20])
ylim([0 lev_deeper])
% plot_rectangle(120,210,-5,5,'b',2)
% plot_rectangle(230,280,-5,5,'r',2)
set(gca,'fontsize',15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_thetaoGlobal_monthly_mean_diff_scaled_d=p_thetaoGlobal_monthly_mean_diff_scaled(find(lev_thetao>=lev_shallower & lev_thetao<=lev_deeper),find(lat_thetao>=lat_lower & lat_thetao<=lat_upper));
for levv=1:numel(lev_1)
    for latt=1:numel(lat_1)
        if p_thetaoGlobal_monthly_mean_diff_scaled_d(levv,latt)<0.1
            plot(lat_1(latt),lev_1(levv),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
        end
    end
end
set(gca,'position',[0.1279    0.1223    0.2377    0.4203])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rectangle(lat_lower_region,lat_upper_region,lev_shallower_region,lev_deeper_region,'k',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.2 Draw scatter
subplot(2,2,4)
h=nan(numel(models_sst)+1,1);
for num=1:numel(models_sst)
    [mar,fco,si,edge]=ew_CMIP_markers_wang(num);
%     if isnan(thetaoGlobal_monthly_mean_diff_region(num))==1
%         h(num)=plot(nan,nan,'marker',mar,'markerfacecolor',fco,'markersize',si+6,'markeredgecolor',edge,'linestyle','none');
%     else
        h(num)=plot(tos_nino34_djf_detrend_std_diff(num),thetaoGlobal_monthly_mean_diff_region(num),'marker',mar,'markerfacecolor',fco,'markersize',si+6,'markeredgecolor',edge,'linestyle','none');
%     end
    hold on
end
h(end)=plot(mean(tos_nino34_djf_detrend_std_diff),mean(thetaoGlobal_monthly_mean_diff_region),'marker','p','markerfacecolor','k','markeredgecolor','k','markersize',25,'linestyle','none');
models_all=cell(numel(models_sst)+1,1);
models_all(1:end-1)=models_sst;
models_all(end)   ={'Multi-model ensemble'};

xlabel('{\Delta}Nino3.4 variability(^oC per degree C of global warming)')
ylabel('SO warming(^oC per degree C of global warming)')
xlim([-0.2 0.4])
ylim([0.2 0.65])
set(gca,'xtick',[-0.2:0.1:0.4],'xticklabel',[-0.2:0.1:0.4])
set(gca,'ytick',[0.2:0.1:0.6],'yticklabel',[0.2:0.1:0.6])
text(-0.2,0.70,'{\bfd} {\Delta}ENSO amplitude {\itvs} SO warming','fontsize',15)
axis normal
set(gca,'tickdir','out')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r rms slope interc n] = ew_regress(tos_nino34_djf_detrend_std_diff,thetaoGlobal_monthly_mean_diff_region);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(tos_nino34_djf_detrend_std_diff,tos_nino34_djf_detrend_std_diff.*slope+interc,'k','linewidth',2)
x=nan(numel(tos_nino34_djf_detrend_std_diff)+2,1);
y=nan(numel(thetaoGlobal_monthly_mean_diff_region)+2,1);
x(2:end-1)=tos_nino34_djf_detrend_std_diff;
y(2:end-1)=thetaoGlobal_monthly_mean_diff_region;
x(1)      =-0.2;
y(1)      =-0.2*slope+interc;
y(end)    =0.2;
x(end)    =(0.2-interc)/slope;
plot(x,x.*slope+interc,'k','linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%
text(-0.19,0.3,['Corre. Coeff.= ' num2str(r,2)],'fontsize',15)
text(-0.19,0.26,['Slope= ' num2str(slope,2)],'fontsize',15)
[r p]=my_corr(tos_nino34_djf_detrend_std_diff(:),thetaoGlobal_monthly_mean_diff_region(:));
if p<0.001
    text(-0.19,0.22,'p<0.001','fontsize',15)
else
    text(-0.19,0.22,['p= ' num2str(p,2)],'fontsize',15)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'fontsize',15)
set(gca,'xminortick','on','yminortick','on')
plot(-100:1:100,ones(numel(-100:1:100),1).*0,'--k')
plot(ones(numel(-100:1:100),1).*0,-100:1:100,'--k')
set(gca,'position',[0.4600    0.1223    0.2377    0.4203])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=tos_nino34_djf_detrend_std_diff;
% y=thetaoGlobal_monthly_mean_diff_region;
% 
% out_vec=fit_ellipse_new(x,y);
% 
% z=[out_vec(1) out_vec(2)];
% plotellipse(z, out_vec(3), out_vec(4), out_vec(5))
legend(h,models_all)











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.3 Correlation between in Qnet and SOI, projected changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_Qnet_onto_SOI=nan(numel(lat_base),numel(lon_base));
% cor_Qnet_onto_SOI=nan(numel(lat_base),numel(lon_base));
% p_Qnet_onto_SOI  =nan(numel(lat_base),numel(lon_base));
% 
% for latt=1:numel(lat_base)
%     for lonn=1:numel(lon_base)
%         d=squeeze(Qnet_monthly_mean_diff(:,latt,lonn));
%         if numel(find(isnan(d)))>=numel(d)*0.5
%         else
%             [r rms slope interc n] = ew_regress(thetaoGlobal_monthly_mean_diff_region,d);
%             reg_Qnet_onto_SOI(latt,lonn)=slope;
%             x=thetaoGlobal_monthly_mean_diff_region;
%             x(find(isnan(d)))=[];
%             d(find(isnan(d)))=[];
%             [r p]=my_corr(x(:),d(:));
%             cor_Qnet_onto_SOI(latt,lonn)=r;
%             p_Qnet_onto_SOI(latt,lonn)  =p;
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat_lower=-80;
% lat_upper=-20;
% lon_left=0;
% lon_right=360;
% 
% lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));
% lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));
% 
% subplot(2,2,2)
% cor_Qnet_onto_SOI_d=cor_Qnet_onto_SOI(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% [c h]=contourf(lon_1,lat_1,cor_Qnet_onto_SOI_d,[-1:0.001:1.0]);
% set(h,'color','none')
% caxis([-1.05 1.05])
% colorbar('Ticks',[-0.9 -0.6 -0.3 0 0.3 0.6 0.9])
% colormap(gca,mm_jet_normal)
% set(gca,'tickdir','out')
% set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
% text(lon_left,lat_upper+5,'{\bfc} Corre. betw. SOI & {\Delta}Qnet, both-scaled','fontsize',15)
% hold on
% coast('k')
% xlim([lon_left lon_right])
% ylim([lat_lower lat_upper])
% set(gca,'fontsize',15)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_Qnet_onto_SOI_d=p_Qnet_onto_SOI(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% for latt=1:numel(lat_1)
%     for lonn=1:numel(lon_1)
%         if p_Qnet_onto_SOI_d(latt,lonn)<0.1
%             plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lon_left_region=0;
% lon_right_region=360;
% lat_lower_region=-70;
% lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.3 Projected changes in zonal wind stress, monthly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('Data_for_Draw_1_v6_SOI_vs_Nino34_Qnet_Tauu_ttest.mat')==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_tauu_monthly_mean_diff_scaled=nan(numel(lat_base),numel(lon_base));
    p_tauu_monthly_mean_diff_scaled=nan(numel(lat_base),numel(lon_base));
    for latt=1:numel(lat_base)
        for lonn=1:numel(lon_base)
            d=squeeze(tauu_monthly_mean_diff(:,latt,lonn));
            if numel(find(isnan(d)))>0
            else
                [h p]=ttest(d,0,'alpha',0.1,'tail','both');
                h_tauu_monthly_mean_diff_scaled(latt,lonn)=h;
                p_tauu_monthly_mean_diff_scaled(latt,lonn)=p;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_tauu_djf_mean_diff_scaled=nan(numel(lat_base),numel(lon_base));
    p_tauu_djf_mean_diff_scaled=nan(numel(lat_base),numel(lon_base));
    for latt=1:numel(lat_base)
        for lonn=1:numel(lon_base)
            d=squeeze(tauu_djf_mean_diff(:,latt,lonn));
            if numel(find(isnan(d)))>0
            else
                [h p]=ttest(d,0,'alpha',0.1,'tail','both');
                h_tauu_djf_mean_diff_scaled(latt,lonn)=h;
                p_tauu_djf_mean_diff_scaled(latt,lonn)=p;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save Data_for_Draw_1_v6_SOI_vs_Nino34_Qnet_Tauu_ttest.mat lon_base lat_base ...
         h_tauu_monthly_mean_diff_scaled p_tauu_monthly_mean_diff_scaled ...
         h_tauu_djf_mean_diff_scaled p_tauu_djf_mean_diff_scaled
else
    load Data_for_Draw_1_v6_SOI_vs_Nino34_Qnet_Tauu_ttest.mat
end


lat_lower=-80;
lat_upper=-20;
lon_left=0;
lon_right=360;

lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));
lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));

subplot(2,2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauu_monthly_mean_diff_mul=nan(numel(lat_base),numel(lon_base));
for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_monthly_mean_diff(:,latt,lonn));
        tauu_monthly_mean_diff_mul(latt,lonn)=mean_withoutnan(d);
    end
end
tauu_monthly_mean_diff_mul_d=tauu_monthly_mean_diff_mul(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
[c h]=contourf(lon_1,lat_1,tauu_monthly_mean_diff_mul_d,[-0.1:0.0001:0.1]);
set(h,'color','none')
caxis([-0.014 0.014])
colorbar('Ticks',[-0.012 -0.008 -0.004 0 0.004 0.008 0.012])
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
%set(gca,'ydir','reverse')
set(gca,'xtick',[0:30:360],'xticklabel',{'','','60^oE','','','','180^oW','','','','60^oW','',''})
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
text(lon_left,lat_upper+5,'{\bfa} Multi-model ensemble mean zonal wind changes','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_tauu_monthly_mean_diff_scaled_d=p_tauu_monthly_mean_diff_scaled(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
for latt=1:3:numel(lat_1)
    for lonn=1:3:numel(lon_1)
        if p_tauu_monthly_mean_diff_scaled_d(latt,lonn)<0.1
            plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
        end
    end
end
set(gca,'position',[0.1285    0.7001    0.2343    0.2653])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon_left_region=0;
lon_right_region=360;
lat_lower_region=-70;
lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8.4 Correlation between in tauu and SOI, projected changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_tauu_onto_SOI=nan(numel(lat_base),numel(lon_base));
cor_tauu_onto_SOI=nan(numel(lat_base),numel(lon_base));
p_tauu_onto_SOI  =nan(numel(lat_base),numel(lon_base));

for latt=1:numel(lat_base)
    for lonn=1:numel(lon_base)
        d=squeeze(tauu_monthly_mean_diff(:,latt,lonn));
        if numel(find(isnan(d)))>=numel(d)*0.5
        else
            [r rms slope interc n] = ew_regress(thetaoGlobal_monthly_mean_diff_region,d);
            reg_tauu_onto_SOI(latt,lonn)=slope;
            x=thetaoGlobal_monthly_mean_diff_region;
            x(find(isnan(d)))=[];
            d(find(isnan(d)))=[];
            [r p]=my_corr(x(:),d(:));
            cor_tauu_onto_SOI(latt,lonn)=r;
            p_tauu_onto_SOI(latt,lonn)  =p;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
% cor_tauu_onto_SOI_d=cor_tauu_onto_SOI(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% [c h]=contourf(lon_1,lat_1,cor_tauu_onto_SOI_d,[-1:0.001:1.0]);
% set(h,'color','none')
% caxis([-1.05 1.05])
% colorbar('Ticks',[-0.9 -0.6 -0.3 0 0.3 0.6 0.9])
reg_tauu_onto_SOI_d=reg_tauu_onto_SOI(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
[c h]=contourf(lon_1,lat_1,reg_tauu_onto_SOI_d,[-1:0.001:1.0]);
set(h,'color','none')
caxis([-0.049 0.049])
colorbar('Ticks',[-0.042 -0.028 -0.014 0 0.014 0.028 0.042])
colormap(gca,mm_jet_normal)
set(gca,'tickdir','out')
set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
set(gca,'xtick',[0:30:360],'xticklabel',{'','','60^oE','','','','180^oW','','','','60^oW','',''})
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
text(lon_left,lat_upper+5,'{\bfc} Zonal wind changes and SO warming','fontsize',15)
hold on
coast('k')
xlim([lon_left lon_right])
ylim([lat_lower lat_upper])
set(gca,'fontsize',15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_tauu_onto_SOI_d=p_tauu_onto_SOI(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
for latt=1:3:numel(lat_1)
    for lonn=1:3:numel(lon_1)
        if p_tauu_onto_SOI_d(latt,lonn)<0.1
            plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon_left_region=0;
lon_right_region=360;
lat_lower_region=-70;
lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)
set(gca,'position',[0.4600    0.7001    0.2343    0.2653])




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 9.1 Correlation between in Qnet and Nino3.4, projected changes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_Qnet_onto_Nino34=nan(numel(lat_base),numel(lon_base));
% cor_Qnet_onto_Nino34=nan(numel(lat_base),numel(lon_base));
% p_Qnet_onto_Nino34  =nan(numel(lat_base),numel(lon_base));
% 
% for latt=1:numel(lat_base)
%     for lonn=1:numel(lon_base)
%         d=squeeze(Qnet_monthly_mean_diff(:,latt,lonn));
%         if numel(find(isnan(d)))>=numel(d)*0.5
%         else
%             [r rms slope interc n] = ew_regress(tos_nino34_djf_detrend_std_diff,d);
%             reg_Qnet_onto_Nino34(latt,lonn)=slope;
%             x=tos_nino34_djf_detrend_std_diff;
%             x(find(isnan(d)))=[];
%             d(find(isnan(d)))=[];
%             [r p]=my_corr(x(:),d(:));
%             cor_Qnet_onto_Nino34(latt,lonn)=r;
%             p_Qnet_onto_Nino34(latt,lonn)  =p;
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat_lower=-80;
% lat_upper=-20;
% lon_left=0;
% lon_right=360;
% 
% lon_1=lon_base(find(lon_base>=lon_left & lon_base<=lon_right));
% lat_1=lat_base(find(lat_base>=lat_lower & lat_base<=lat_upper));
% 
% figure
% subplot(2,2,1)
% cor_Qnet_onto_Nino34_d=cor_Qnet_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% [c h]=contourf(lon_1,lat_1,cor_Qnet_onto_Nino34_d,[-1:0.001:1.0]);
% set(h,'color','none')
% caxis([-1.05 1.05])
% colorbar('Ticks',[-0.9 -0.6 -0.3 0 0.3 0.6 0.9])
% colormap(gca,mm_jet_normal)
% set(gca,'tickdir','out')
% set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
% text(lon_left,lat_upper+5,'{\bfa} Corre. betw. Nino34 & {\Delta}Qnet, both-scaled','fontsize',15)
% hold on
% coast('k')
% xlim([lon_left lon_right])
% ylim([lat_lower lat_upper])
% set(gca,'fontsize',15)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_Qnet_onto_Nino34_d=p_Qnet_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% for latt=1:numel(lat_1)
%     for lonn=1:numel(lon_1)
%         if p_Qnet_onto_Nino34_d(latt,lonn)<0.1
%             plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lon_left_region=0;
% lon_right_region=360;
% lat_lower_region=-70;
% lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 9.1 Correlation between in Qnet and Nino3.4, projected changes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_tauu_onto_Nino34=nan(numel(lat_base),numel(lon_base));
% cor_tauu_onto_Nino34=nan(numel(lat_base),numel(lon_base));
% p_tauu_onto_Nino34  =nan(numel(lat_base),numel(lon_base));
% 
% for latt=1:numel(lat_base)
%     for lonn=1:numel(lon_base)
%         d=squeeze(tauu_monthly_mean_diff(:,latt,lonn));
%         if numel(find(isnan(d)))>=numel(d)*0.5
%         else
%             [r rms slope interc n] = ew_regress(tos_nino34_djf_detrend_std_diff,d);
%             reg_tauu_onto_Nino34(latt,lonn)=slope;
%             x=tos_nino34_djf_detrend_std_diff;
%             x(find(isnan(d)))=[];
%             d(find(isnan(d)))=[];
%             [r p]=my_corr(x(:),d(:));
%             cor_tauu_onto_Nino34(latt,lonn)=r;
%             p_tauu_onto_Nino34(latt,lonn)  =p;
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,2,3)
% cor_tauu_onto_Nino34_d=cor_tauu_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% [c h]=contourf(lon_1,lat_1,cor_tauu_onto_Nino34_d,[-1:0.001:1.0]);
% set(h,'color','none')
% caxis([-1.05 1.05])
% colorbar('Ticks',[-0.9 -0.6 -0.3 0 0.3 0.6 0.9])
% colormap(gca,mm_jet_normal)
% set(gca,'tickdir','out')
% set(gca,'ytick',[-80:10:80],'yticklabel',{'80^oS','','60^oS','','40^oS','','20^oS','','0^o','','20^oN','','40^oN','','60^oN','','80^oN'})
% set(gca,'xtick',[0:30:360],'xticklabel',{'0^o','','60^oE','','120^oE','','180^oW','','120^oW','','60^oW','',''})
% text(lon_left,lat_upper+5,'{\bfa} Corre. betw. Nino34 & {\Delta}tauu, both-scaled','fontsize',15)
% hold on
% coast('k')
% xlim([lon_left lon_right])
% ylim([lat_lower lat_upper])
% set(gca,'fontsize',15)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_tauu_onto_Nino34_d=p_tauu_onto_Nino34(find(lat_base>=lat_lower & lat_base<=lat_upper),find(lon_base>=lon_left & lon_base<=lon_right));
% for latt=1:numel(lat_1)
%     for lonn=1:numel(lon_1)
%         if p_tauu_onto_Nino34_d(latt,lonn)<0.1
%             plot(lon_1(lonn),lat_1(latt),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',1)
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lon_left_region=0;
% lon_right_region=360;
% lat_lower_region=-70;
% lat_upper_region=-50;
% plot_rectangle(lon_left_region,lon_right_region,lat_lower_region,lat_upper_region,'k',2)
% 


































