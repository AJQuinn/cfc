function [dat] = normdata(dat)

    dat = (dat - min(dat)) / ( max(dat) - min(dat) );
