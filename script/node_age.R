# Node age range among environmental types

node.age.df <- gdata$data[pmatch(gdata$data$label[1:2928],V(conet)$name),]
node.age.df1 <- rbind(data.frame(node.age.df[V(conet)$s1 == 1,c(6,3,8)],env=env.name1[1]),
                      data.frame(node.age.df[V(conet)$s2 == 1,c(6,3,8)],env=env.name1[2]),
                      data.frame(node.age.df[V(conet)$s3 == 1,c(6,3,8)],env=env.name1[3]),
                      data.frame(node.age.df[V(conet)$s4 == 1,c(6,3,8)],env=env.name1[4]),
                      data.frame(node.age.df[V(conet)$s5 == 1,c(6,3,8)],env=env.name1[5]),
                      data.frame(node.age.df[V(conet)$s6 == 1,c(6,3,8)],env=env.name1[6]),
                      data.frame(node.age.df[V(conet)$s7 == 1,c(6,3,8)],env=env.name1[7]),
                      data.frame(node.age.df[V(conet)$s8 == 1,c(6,3,8)],env=env.name1[8]),
                      data.frame(node.age.df[V(conet)$s9 == 1,c(6,3,8)],env=env.name1[9]),
                      data.frame(node.age.df[V(conet)$s10 == 1,c(6,3,8)],env=env.name1[10]),
                      data.frame(node.age.df[V(conet)$s11 == 1,c(6,3,8)],env=env.name1[11]),
                      data.frame(node.age.df[V(conet)$s12 == 1,c(6,3,8)],env=env.name1[12]),
                      data.frame(node.age.df[V(conet)$s13 == 1,c(6,3,8)],env=env.name1[14]),
                      data.frame(node.age.df[V(conet)$s14 == 1,c(6,3,8)],env=env.name1[13])
                      )

branch.mean <- node.age.df1 %>%
        mutate(age = branch.length + branch) %>%
                  group_by(env) %>%
                  summarise(b.mean = quantile(log10(age),probs = .5)) %>%
                  arrange(desc(b.mean))
        
node.age.df1$env <- factor(node.age.df1$env,levels = branch.mean$env[c(1,2,8:12,14,3:7,13)])

ggplot(node.age.df1, aes(x=env,y=branch.length,color=env))+
        geom_quasirandom(width = .3,size=1)+
        scale_y_log10(breaks = c(0.001,0.01,.1),
                      label=expression(10^-3,10^-2,10^-1))+
        geom_vline(xintercept = 8.5, linetype = 2)+
        theme_bw()+guides(color=FALSE)+
        scale_color_discrete_qualitative(alpha = .4)+
        geom_violin(fill=NA,draw_quantiles = 0.5,color=grey(.2))+
        theme(axis.text.x = element_text(angle=-75,hjust = 0))+
        xlab("Environment type")+ylab("Scaled vertex age")
ggsave("age_env.pdf",height = 5,width = 5)
