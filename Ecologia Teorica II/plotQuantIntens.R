#Function to plot partial responses

plotQuantIntens <- function(imList, noCuts = NULL, Quad, p.pp, dir = "", name = ""){

    require(spatstat)
    
    Z <- imList
    
    if(is.null(noCuts)){
        pdf(paste0(dir, name, ".pdf"), width = 5, height = 5)
        for(i in seq_along(Z)){
            plot(rhohat(Quad$data, Z[[i]]) , main = paste(name, names(Z)[i], sep = ", "))
        }
        dev.off()
    } else {
        cuts <- seq(0, 1, by = 1/noCuts)
        
        quants <- lapply(Z, function(x){seq(min(x), max(x), len = noCuts + 1)})
        
        Z.cut <- foreach(k = seq_along(Z)) %do% {
            cut(Z[[k]], breaks = quants[[k]], labels = 1:(length(cuts)-1))
        }
        
        V <- lapply(Z, function(x)tess(image = x))
        
        counts <- foreach(k = seq_along(V)) %do% {
            quadratcount(Quad$data, tess = V[[k]])
        }
        
        pdf(paste0(dir, name, ".pdf"), width = 15, height = 5)
        par(mfrow = c(1,3))
        for(i in seq_along(Z)){
            plot(counts[[i]], main = paste(name, names(Z)[i], sep = ", "))
            plot(V[[i]], main = "")
            points(p.pp, pch = "+", col = "green", cex = 2)
            plot(rhohat(Quad$data, Z[[i]]) , main = "")
        }
        dev.off()   
    }
}