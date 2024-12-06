ggplot() +
  
  geom_line(data=data1,aes(x=log(data1$V1,10),  y=data1$bf*data1$V1) , lwd=0.5)+
  geom_line(data=data2,aes(x=log(data2$V1,10),  y=data2$bf*data2$V1),color="red" , lwd=0.5)
#red is fungi