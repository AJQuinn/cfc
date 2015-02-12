function [dat] = norm_dat(dat)

    dat = (dat - min(dat)) / ( max(dat) - min(dat) );
