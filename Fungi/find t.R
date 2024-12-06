ggplot(p, aes(x = T))+geom_vline(aes(xintercept=550),linetype="dashed")+
  geom_hline(aes(yintercept=106/139),linetype="dashed",color="cyan")+
  geom_line(aes(y = (¦¤GCF+106)/139, color = "¦¤GCF"))+
  geom_line(aes(y = Vscore, color = "Vscore")) +scale_y_continuous(limits =c(0,1),sec.axis = sec_axis(~.*139-106),name = "")+
  scale_x_continuous(limits = c(100,1500),breaks=c(100,300,500,700,900,1100,1300,1500))+scale_color_manual(values = c("pink","cyan"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 12))+theme_bw()+scale_fill_brewer()+annotate("text",x=550,y=0,label="T=550")+annotate("text",x=550,y=0.71,label="¦¤GCF=-2")+annotate("text",x=550,y=0.9,label="Vscore=0.92")
ggsave("D:\\bigslice\\t.tiff",width = 5,height = 4)
