emp_GASP_plot <-
function(em, fun, data, emul_type, exp.ql, exp.qu,labels, ylab, xlab,ylim,col_CI_area,col_points,col_fun,col_mean, points = FALSE) {
  xax = data$training
  if (missing(points)) {
    points = FALSE
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
  f.values = fun(input)
  coverage = f.values <= qu & f.values >= ql
  if (any(coverage == FALSE)) {
    for (i in which(coverage == FALSE)) {
      coverage[i] = all.equal(f.values[i],ql[i])==TRUE || all.equal(f.values[i],qu[i])==TRUE
    }
  }
  emp.coverage = sum(coverage)/length(f.values)
  if (missing(ylim)) ylim = range(f.values,ql,qu)
  plot(input,ql,type = "l", col = col_CI_area, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
  par(new = TRUE)
  plot(input,qu,type = "l", col = col_CI_area, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
  polygon(c(input,rev(input)),c(qu,rev(ql)),col=col_CI_area,border = NA)
  
  par(new = TRUE)
  plot(input,em$mu,type = "l",lty = 1, col=col_mean, lwd = 3, ylab = "",xlab ="",ylim = ylim, xaxt = "n", yaxt = "n")
  par(new = TRUE)
  plot(input,f.values,type = "l", lty = 1,col = col_fun, lwd = 3, ylab = "", xlab = "",ylim = ylim, xaxt = "n", yaxt = "n")
  par(new = TRUE)
  
  plot(input,exp.ql,type = "l", xaxt = "n",lwd = 1, lty = 4, col = "black", ylim = ylim, ylab = "",xlab = "",xaxt = "n", yaxt = "n")
  par(new = TRUE)
  plot(input,exp.qu,type = "l", xaxt = "n",lwd = 1, lty = 4, col = "black", ylim = ylim, ylab = "",xlab = "",xaxt = "n", yaxt = "n")
  
  if (points) points(xax,fun(xax),pch = 20,col = col_points,cex = 2.5)
  axis(1, at=xax, labels=labels, cex.axis = 1.5)
  axis(2, cex.axis = 1.5)
  
  yl = deparse(ylab)
  
  legend("topright",lty = c(rep(1,3),NA),pch=c(rep(NA,3),20),lwd = rep(3,3),cex = 1.4,
         bty = "n",
         legend = c("True function","Emulator mean","95%-CI area","Observed"),
         col = c(col_fun,col_mean,col_CI_area,col_points))
  title(main = "",ylab = ylab, xlab = xlab,cex.lab = 2, mgp=c(3,1,.5))
  #text(sum(range(input))/2,ylim[1],bquote(paste("Emp. coverage "==.(round(emp.coverage,4)*100),"%")))
  par(new = FALSE)
}
