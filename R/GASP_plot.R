GASP_plot <-
  function(em, fun, data, emul_type, labels, yax, ylab, xlab,ylim,col_CI_area,col_points,col_fun,col_mean,plot_training = FALSE, plot_fun = TRUE) {
    xax = data$training
    if(missing(fun)) {
      plot_fun = FALSE
    }
    if(missing(yax)) {
      yax = c(-1,0,1)
    }
    if(missing(col_CI_area)) {
      col_CI_area = "lightgreen"
    }
    if(missing(col_points)) {
      col_points = "darkorchid"
    }
    if(missing(col_fun)) {
      col_fun = "navajowhite"
    }
    if(missing(col_mean)) {
      col_mean = "darkgreen"
    }
    if(missing(labels)) {
      labels <- xax
    }
    if(missing(xlab)) {
      xlab <- "y, input"
    }
    if(missing(ylab)) {
      ylab <- "g, output"
    }
    input = em$x; var = diag(em$var)
    ql = em$mu + qnorm(.025,0,1)*sqrt(var)
    qu = em$mu + qnorm(.975,0,1)*sqrt(var)
    if (!missing(fun)) f.values = fun(input)
    ylim = range(ql,qu)
    if (exists("f.values")) {ylim = range(f.values,ylim)}
    plot(input,ql,type = "l", col = col_CI_area, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
    par(new = TRUE)
    plot(input,qu,type = "l", col = col_CI_area, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
    polygon(c(input,rev(input)),c(qu,rev(ql)),col=col_CI_area,border = NA)
    
    par(new = TRUE)
    plot(input,em$mu,type = "l",lty = 1, col=col_mean, lwd = 3, ylab = "",xlab ="",ylim = ylim, xaxt = "n", yaxt = "n")
    par(new = TRUE)
    if (plot_fun) {
      plot(input,f.values,type = "l", lty = 1,col = col_fun, lwd = 3, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
      par(new = TRUE)
    }
    
    if (plot_training) {points(data$training,data$fD,pch = 20,col = col_points,cex = 2.5)}
    axis(1, at=xax, labels=labels, cex.axis = 2)
    axis(2, at = yax,cex.axis = 2)
    
    yl = deparse(ylab)
    # legend("topleft", lty = 1, lwd = 3, cex = 1.4,
    #        bty = "n",
    #        legend = c("True function"),
    #        col = c(col_fun))
    
    # legend("topleft",lty = c(rep(1,1),NA),pch=c(rep(NA,1),20),lwd = rep(3,1),cex = 1.4,
    #        bty = "n",
    #        legend = c("True function","Observed"),
    #        col = c(col_fun,col_points))
    
    
    # legend("topright",lty = c(rep(1,3),NA),pch=c(rep(NA,3),20),lwd = rep(3,3),cex = 1.4,
    #        bty = "n",
    #        legend = c("True function","Emulator mean","95%-CI area","Observed"),
    #        col = c(col_fun,col_mean,col_CI_area,col_points))
    title(main = emul_type,ylab = ylab, xlab = xlab,cex.lab = 2.5, mgp=c(2.8,1.5,.5))
    par(new = FALSE)
  }