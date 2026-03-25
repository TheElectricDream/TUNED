function [removed_pixels] = recoverRemovedEvents(filter_mask, counts)


    temp = filter_mask;
    temp(filter_mask>0)=0;
    temp(filter_mask==0)=1;
    counts_temp = counts.*temp;
    removed_pixels = (counts_temp>0).*1.0;
    
end