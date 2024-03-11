rm(list = ls(all=TRUE))
library(ggplot2)
df_cytb = read.csv('Midori2_birds_cytb_ghahskew_better.csv')
df_cytb = df_cytb[,c(2,3,4)]
df_cytb1 <- df_cytb[order(df_cytb$GhAhSkew_seq_begining, decreasing = TRUE), ] 
df_cytb1 = df_cytb1[,c(1,2)]
df_cytb2 <- df_cytb[order(df_cytb$GhAhSkew_start_codon, decreasing = TRUE), ]
df_cytb2 = df_cytb2[,c(1,3)]
df_cytb1_1 <- df_cytb[order(df_cytb$GhAhSkew_seq_begining, decreasing = FALSE), ] 
df_cytb1_1 = df_cytb1_1[,c(1,2)]
df_cytb2_1 <- df_cytb[order(df_cytb$GhAhSkew_start_codon, decreasing = FALSE), ]
df_cytb2_1 = df_cytb2_1[,c(1,3)]
df_cytb1$class = 'Aves'
df_cytb2$class = 'Aves'
df_cytb1_1$class = 'Aves'
df_cytb2_1$class = 'Aves'
df_cytb1_percent_max <- Reduce(rbind,                                 # Top N highest values by group
                    by(df_cytb1,
                       df_cytb1["class"],
                       head,
                       n = 72))
df_cytb2_percent_max <- Reduce(rbind,                                 # Top N highest values by group
                           by(df_cytb2,
                              df_cytb2["class"],
                              head,
                              n = 72))
df_cytb1_1_percent_max <- Reduce(rbind,                                 # Top N highest values by group
                               by(df_cytb1_1,
                                  df_cytb1_1["class"],
                                  head,
                                  n = 72))
df_cytb2_2_percent_max <- Reduce(rbind,                                 # Top N highest values by group
                               by(df_cytb2_1,
                                  df_cytb2_1["class"],
                                  head,
                                  n = 72))

ggplot(df_cytb1_1_percent_max, aes(x = class, y = GhAhSkew_seq_begining))+
  geom_boxplot()



ggp <- ggplot(NULL, aes(class, GhAhSkew_seq_begining)) +    # Draw ggplot2 plot based on two data frames
  geom_boxplot(data = df_cytb1_percent_max, col = "red") +
  geom_boxplot(data = df_cytb1_1_percent_max, col = "blue")
ggp    
ggp1 <- ggplot(NULL, aes(class, GhAhSkew_start_codon)) +    # Draw ggplot2 plot based on two data frames
  geom_boxplot(data = df_cytb1_percent_max, col = "red") +
  geom_boxplot(data = df_cytb1_1_percent_max, col = "blue")
ggp1  
write.csv(df_cytb1_percent_max, file = 'Start_seq_max.csv')
write.csv(df_cytb2_percent_max, file = 'Start_codon_max.csv')
write.csv(df_cytb1_1_percent_max, file = 'Start_seq_min.csv')
write.csv(df_cytb2_2_percent_max, file = 'Start_codon_min.csv')
