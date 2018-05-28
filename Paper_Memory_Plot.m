clear; close all;

for size_i = 1:8
    file_name = sprintf("paper_timing_%.2d",size_i);
    load(file_name);
    temp = whos('Rmat');
    memMVM_array(size_i) = temp.bytes;
    temp = whos('CovKernel');
    memFil_array(size_i) = temp.bytes;
end
   
save paper_memory_all memMVM_array memFil_array
WFS_size = (1:8)*6+5;

plot(WFS_size,memMVM_array);
hold on
plot(WFS_size,memFil_array);