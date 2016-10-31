## figure instead of Table 3 and 4
library(ggplot2)
psi10 <- c(0,0,0,0,0,0,0,0,0)

# this is information taken from Table 4
widthN <- c(0.546, 0.545, 0.590, 0.597, 0.625, 0.588, 0.663, 0.598, 0.566)
width0.05 <- c(0.624, 0.625, 0.655, 0.639, 0.681, 0.633, 0.707, 0.639, 0.641)
width0.1 <- c(0.702, 0.701, 0.721, 0.680, 0.737, 0.680, 0.752, 0.681, 0.718)
widthA <- c(0.618, 0.615, 0.666, 0.637, 0.683, 0.637, 0.713, 0.638, 0.648)

positionN <- seq(1,29,by=3.5)
scenario <- seq(1,9)
dataStart <- as.data.frame(cbind(posN = rep(positionN, each=2), pos0.05=rep(positionN+0.5, each=2), pos0.1=rep(positionN+1, each=2), posA=rep(positionN+1.5, each=2), sc = rep(scenario, each=2), extN=c(rbind(psi10-widthN/2,psi10+widthN/2)), ext0.05=c(rbind(psi10-width0.05/2,psi10+width0.05/2)), ext0.1=c(rbind(psi10-width0.1/2,psi10+width0.1/2)), extA=c(rbind(psi10-widthA/2,psi10+widthA/2))))

id <- seq(1,72*2)
yNup <- rep(c(positionN, positionN+0.5, positionN+1, positionN+1.5)+0.15,each=2)
yNlo <- rep(c(positionN, positionN+0.5, positionN+1, positionN+1.5)-0.15,each=2)
xNlo <- c(psi10-widthN/2, psi10-width0.05/2, psi10-width0.1/2, psi10-widthA/2)
xNup <- c(psi10+widthN/2, psi10+width0.05/2, psi10+width0.1/2, psi10+widthA/2)
# this is information taken from Table 4
cov <- c(0.926, 0.944, 0.923, 0.935, 0.941, 0.935, 0.921, 0.935, 0.889,0.970, 0.967, 0.958, 0.954, 0.939, 0.949, 0.954, 0.956, 0.939, 0.985, 0.970, 0.968, 0.965, 0.963, 0.949, 0.958, 0.969, 0.974, 0.965, 0.970, 0.953, 0.955, 0.954, 0.953 ,0.942, 0.958, 0.926)
id_italic <- c(1,3,4,6,7,8,9,10,11,19,20,21,22,26,27,28,29,36)
tickMar <- as.data.frame(cbind(gr=rep(id,each=2), posX=rep(c(rbind(xNlo,xNup)),each=2), posY=c(rbind(yNup,yNlo))))

ggplot(data=dataStart, aes(x=extN,y=posN,group=sc))+
  geom_line(aes(fill="nn"), size=0.5, show_guide=TRUE)+
  geom_line(aes(x=ext0.05,y=pos0.05,group=sc, fill="mn 0.05"), size=0.5, show_guide=TRUE, linetype="dashed")+
  geom_line(aes(x=ext0.1,y=pos0.1,group=sc, fill="mn 0.1"), size=0.5, show_guide=TRUE, linetype="dotted")+
  geom_line(aes(x=extA,y=posA,group=sc, fill="mn adaptive"), size=0.5, show_guide=TRUE, linetype="twodash")+
  geom_line(data=tickMar, aes(x=posX,y=posY,group=gr), size=0.3)+
  #xlab(expression(psi[10]))+
  ylab("")+
  scale_x_continuous(name=expression(psi[10]), limits=c(-0.9,0.9),breaks=0, labels="")+
  scale_y_continuous(breaks=seq(0.5,32.5,by=3.5), labels=c("","Sc.1 NR", "Sc.2 NNR", "Sc.3 NR", "Sc.4 NNR", "Sc.5 NR", "Sc.6 R", "Sc.7 R", "Sc.8 NR", "Sc.9 NNR"))+
  theme(axis.text.y=element_text(angle=0, vjust=3.8), legend.text=element_text())+
  scale_fill_manual(
    "CI method", values=rep(1,4),
    guide=guide_legend(override.aes = list(linetype=c("dashed", "dotted", "twodash", "solid")))
  )+
  annotate("text", x=xNup[id_italic]+0.07,y=c(positionN,positionN+0.5,positionN+1,positionN+1.5)[id_italic],label=cov[id_italic], size=2.5, fontface=2)+
  annotate("text", x=xNup[-id_italic]+0.07,y=c(positionN,positionN+0.5,positionN+1,positionN+1.5)[-id_italic],label=cov[-id_italic], size=2.5)


