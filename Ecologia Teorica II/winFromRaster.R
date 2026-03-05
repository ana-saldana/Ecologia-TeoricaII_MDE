winFromRaster <- function(r){      
    require(spatstat); require(foreach)
    
    r.df <- as.data.frame(r, xy = T)
    
    ux = sort(unique(r.df$x)) #Extracting unique coordinates
    uy = sort(unique(r.df$y))
    nx = length(ux) #length of unique coordinates
    ny = length(uy)
    ref.cols = match(r.df$x, ux) #position of every data point
    ref.lines = match(r.df$y, uy)
    vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
    ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
    vec[ref.vec] = 1
    data.mask = matrix(vec, max(ref.lines), max(ref.cols))
    w = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window
    
    return(w)
}
