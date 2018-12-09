TGASP_plot <-
function(tem,fun,data,labels,ylim,points = FALSE) {
  xax = data$training
  if(missing(labels)) {
    labels <- xax
  }
  input = tem$x; sigma = diag(tem$sigma); df = tem$df
  ql = tem$mu + qnorm(.025,0,1)*sqrt(sigma)
  qu = tem$mu + qnorm(.975,0,1)*sqrt(sigma)
  f.values = fun(input)
  # coverage = f.values <= qu & f.values >= ql
  # if (any(coverage == FALSE)) {
  #   for (i in which(coverage == FALSE)) {
  #     coverage[i] = all.equal(f.values[i],ql[i])==TRUE || all.equal(f.values[i],qu[i])==TRUE
  #   }
  # }
  # emp.coverage = sum(coverage)/length(f.values)
  if (missing(ylim)) ylim = range(f.values,ql,qu)
  plot(input,f.values,type = "l", lty = 1,col = "black", ylab = "", xlab = "",ylim = ylim, xaxt = "n")
  par(new = TRUE)
  plot(input,tem$mu,type = "l",lty = 5,col="black",ylab = "",xlab ="",ylim = ylim, xaxt = "n")
  par(new = TRUE)
  plot(input,ql,type = "l", col = "black", lty = 3,ylab = "", xlab = "",ylim = ylim, xaxt = "n")
  par(new = TRUE)
  plot(input,qu,type = "l", col = "black", lty = 3,ylab = "", xlab = "",ylim = ylim, xaxt = "n")
  axis(1, at=xax, labels=labels)
  points(xax,fun(xax),pch = 1)
  legend("topright",lty=c(1,5,3),
         bty = "n",
         legend = c("True function","Emulator mean","95%-CIs"),col = c("black","black","black"))
  title(main = "T-GASP",ylab = "f(x)", xlab = "x",cex.lab = 1.3)
  # text(sum(range(input))/2,ylim[1],bquote(paste("Emp. coverage "==.(round(emp.coverage,2)*100),"%")))
  par(new = FALSE)
}
