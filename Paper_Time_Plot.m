clear; close all;

for size_i = 1:8
    file_name = sprintf("paper_timing_%.2d",size_i);
    load(file_name);
    timeMVM_array(size_i) = Time_MVM;
    timeFil_array(size_i) = Time_Fil;
end
   
save paper_timing_all timeMVM_array timeFil_array
WFS_size = (1:8)*6+5;

plot(WFS_size,timeMVM_array);
hold on
plot(WFS_size,timeFil_array);