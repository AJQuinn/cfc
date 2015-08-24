function [dat] = cfc_util_normdata(dat)

    dat = (dat - min(dat)) / ( max(dat) - min(dat) );
