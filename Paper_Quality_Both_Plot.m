clear; close all;

sigma_2 = 0.148;
size_i = 4;
WFS_size = size_i*6+5

file_name = sprintf("paper_quality_noise%.3f_size%.2d.mat",sigma_2,size_i);
load(file_name);


%****************************************************
%WFE
Time_MVM_avg = mean(Time_MVM);
WFE_layers_MVM_avg = mean(WFE_layers_MVM);
WFE_Target_MVM_avg = mean(WFE_Target_MVM);
WFE_layers_MVM_var = var(WFE_layers_MVM);
WFE_Target_MVM_var = var(WFE_Target_MVM);

ITER = mean([ITER_X ITER_Y],2);
Time_Fil_avg = mean(Time_Fil,2);
WFE_layers_FilterNI_avg = mean(WFE_layers_FilterNI,2);
WFE_Target_FilterNI_avg = mean(WFE_Target_FilterNI,2);
WFE_layers_FilterNI_var = var(WFE_layers_FilterNI,[],2);
WFE_Target_FilterNI_var = var(WFE_Target_FilterNI,[],2);


plot(Time_Fil_avg,WFE_layers_FilterNI_avg,'Marker','*','LineStyle','none');
hold on
plot(Time_MVM_avg,WFE_layers_MVM_avg,'Marker','o','LineStyle','none');
text(Time_Fil_avg,WFE_layers_FilterNI_avg,num2str(round(ITER)))
ylabel('WFE (nm)');
xlabel('Time (s)');

figure;
plot(Time_Fil_avg,WFE_Target_FilterNI_avg,'Marker','*','LineStyle','none');
hold on
plot(Time_MVM_avg,WFE_Target_MVM_avg,'Marker','o','LineStyle','none');
text(Time_Fil_avg,WFE_Target_FilterNI_avg,num2str(round(ITER)))
ylabel('WFE (nm)');
xlabel('Time (s)');

%****************************************************
% X Slope Error

Error_TargetX_FilterNI_avg = mean(Error_TargetX_FilterNI,2);
Error_layersX_FilterNI_avg = mean(Error_layersX_FilterNI,2);
Error_TargetX_MVM_avg = mean(Error_TargetX_MVM);
Error_layersX_MVM_avg = mean(Error_layersX_MVM);

figure;
plot(Time_Fil_avg,Error_layersX_FilterNI_avg,'Marker','*','LineStyle','none');
hold on
plot(Time_MVM_avg,Error_layersX_MVM_avg,'Marker','o','LineStyle','none');
text(Time_Fil_avg,Error_layersX_FilterNI_avg,num2str(round(ITER)))
ylabel('Slope Error');
xlabel('Time (s)');

figure;
plot(Time_Fil_avg,Error_TargetX_FilterNI_avg,'Marker','*','LineStyle','none');
hold on
plot(Time_MVM_avg,Error_TargetX_MVM_avg,'Marker','o','LineStyle','none');
text(Time_Fil_avg,Error_TargetX_FilterNI_avg,num2str(round(ITER)))
ylabel('Slope Error');
xlabel('Time (s)');