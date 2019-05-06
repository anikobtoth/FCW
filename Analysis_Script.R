library(tidyverse)
library(vegan)
library(reshape2)
library(cowplot)

source('~/Desktop/MacQuarie_PhD/Thesis/Chapter1_SimPairs/FCW/Pair_Functions.R')

load("EastAfricanMammalData/speciesData.RData")
load("EastAfricanMammalData/siteData.RData")
load("EastAfricanMammalData/Occurrence_data_longformat.RData")

### PA tables ####
PA <- dcast(BINOMIAL~sitekey, data = data, value.var = 'observed')
PA <- namerows(PA)
pa <- PA[which(rownames(PA) %in% sppdat$BINOMIAL[sppdat$MASS_KG >= 1]),]
PA <- tobinary(list(PA, pa)) %>% map(clean.empty)

pa <- PA[[2]] # species > 1 kg
PA <- PA[[1]] # all species

### species split by broad diet ####
#out <- simpairs(PA) %>% dist2edgelist(PA)
out <- FETmP_Pairwise(pa) %>% dist2edgelist(pa) #slow step
out <- out[out$Sp1 != out$Sp2,]
out$diet.Sp1 <- sppdat[out$Sp1,"DIET1.5"]
out$diet.Sp2 <- sppdat[out$Sp2,"DIET1.5"]

out$mass.Sp1 <- sppdat[out$Sp1,"MASS_KG"]
out$mass.Sp2 <- sppdat[out$Sp2,"MASS_KG"]

out$Z.Score <- qnorm(out$Score)

out$diet.pair <- map2(out$diet.Sp1, out$diet.Sp2, function(x, y) c(x,y)) %>% map(sort) %>% map(paste, collapse = "-") %>% unlist()
out$diet.match <- out$diet.Sp1 == out$diet.Sp2

out$type <- posnegzero(out$Z.Score)
out %>% group_by(diet.pair) %>% summarise(agg = percpos(Z.Score))

out2 <- out[!out$type == "ZERO"  & out$diet.match == TRUE,]
ggplot(out2, aes(y = pnorm(Z.Score), x = diet.pair, fill = diet.pair)) + 
  geom_boxplot(notch = T) + facet_grid(type~., scales = "free") + 
  scale_fill_hue(h = c(15, 215, 375), l = c(50, 80, 60), c = c(50, 90, 80)) + 
  theme(axis.text = element_text(size = 12), axis.title.x = element_blank(), 
        legend.position = "none", panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  scale_x_discrete(labels = c("carnivores", "herbivores", "omnivores")) + labs(y = "Co-occurrence probability") + 
  panel_border(remove=F, colour = "black")

out2 %>% group_by(diet.pair) %>% summarise(agg = percpos(Z.Score), count = length(Z.Score))


### pair split####

out <- FETmP_Pairwise(pa) %>% dist2edgelist(pa)
out <- out[out$Sp1 != out$Sp2,]
out$Z.Score <- qnorm(out$Score)
out$diet.Sp1 <- sppdat[out$Sp1,"DIET2"]
out$diet.Sp2 <- sppdat[out$Sp2,"DIET2"]

out$diet.match <- out$diet.Sp1 == out$diet.Sp2
out$type <- posnegzero(out$Z.Score)
out$score <- pnorm(out$Z.Score)*2-1

ggplot(out[!out$type == "ZERO", ], aes(y = pnorm(Z.Score), x = diet.match, fill = diet.match)) + 
  geom_boxplot(notch = T) + facet_grid(type~., scales = "free") + 
  scale_x_discrete(labels = c("Different", "Same")) +
  theme(axis.text = element_text(size = 12), legend.position = "none", 
        panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  scale_fill_hue(h = c(90, 270), l = c(75, 50)) + labs(y = "Co-occurrence probability", x = "Pair diet") + panel_border(remove = F, col = "black")

out %>% group_by(diet.match) %>% summarise(agg = percpos(Z.Score), count = length(Z.Score))
out %>% group_by(diet.match, type) %>% summarise(median = median(pnorm(Z.Score)))

##### Time-based site split in pa #####

pas <- list(pa[,which(colnames(pa) %in% sitedat$sitekey[sitedat$timeybp <= 50])], 
            pa[,which(colnames(pa) %in% sitedat$sitekey[sitedat$timeybp > 50])])

pas <- map(pas, clean.empty, minrow = 1)
names(pas) <- c("Recent", "Historical")

out <- pas %>% map(FETmP_Pairwise) %>% map2(pas, dist2edgelist) %>% bind_rows(.id = "time")
out <- out[out$Sp1 != out$Sp2,]
out$time <- factor(out$time, levels = c("Recent", "Historical"))
out$Z.Score <- qnorm(out$Score)
out$type <- posnegzero(out$Z.Score)

out$diet.Sp1 <- sppdat[out$Sp1,"DIET1"]
out$diet.Sp2 <- sppdat[out$Sp2,"DIET1"]

out$diet.pair <- map2(out$diet.Sp1, out$diet.Sp2, function(x, y) c(x,y)) %>% map(sort) %>% map(paste, collapse = "-") %>% unlist()
out$diet.match <- out$diet.Sp1 == out$diet.Sp2

## Overall difference between time intervals
b <- ggplot(out[out$type != "ZERO",], aes(y = pnorm(Z.Score), x = time, fill = time)) + 
  geom_boxplot(notch = T) + facet_grid(type~., scales = "free") + labs(x = "Time", y = "Co-occurrence probability") +
  scale_fill_hue(h = c(45, 250), l = c(80, 50), c = c(100, 10)) + 
  theme(legend.position = "none", panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        plot.margin = unit(c(0.5, .5, 0, 0), "cm")) +panel_border(remove = FALSE, col = 'black')

a <- ggplot(out, aes(x = pnorm(Z.Score), col = time)) + geom_density() + coord_flip() +
  scale_colour_hue(h = c(45, 250), l = c(80, 50), c = c(100, 10)) + 
  theme(legend.position = 'none',panel.background = element_blank(), panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12), plot.margin = unit(c(0.5, 0, 0.5, .5), "cm")) + 
  labs(x = "Co-occurrence probability", y = "Relative frequency") 

plot_grid(a+panel_border(remove = FALSE, col = 'black'), 
          b + theme(axis.title.y = element_blank(), panel.spacing = unit(0, "lines")), 
          align = "hv", rel_widths = c(3,5), axis = "tb")

# Time interval side by side
ggplot(out[out$type != "ZERO",], aes(y = pnorm(Z.Score), x = time, fill = time)) + 
  geom_boxplot(notch = T) + facet_grid(type~diet.match, scales = 'free') +
  scale_fill_hue(h = c(45, 250), l = c(80, 50), c = c(100, 10)) + labs(y = "Co-occurrence probability", x = "Time interval") +
  theme(axis.text = element_text(size = 12), legend.position = "none", 
        panel.background = element_blank(), panel.border = element_rect(fill = NA)) + panel_border(remove = F, col = "black")


# Fine diet categories
ggplot(out[out$type != "ZERO",], aes(y = pnorm(Z.Score), x = diet.pair, fill = diet.pair)) + 
  geom_boxplot(notch = T) + facet_grid(type~time, scales = 'free') + scale_x_discrete(labels = c("a_dom", "a-p_dom", "p_dom")) +
  scale_fill_hue(h = c(120, 250, 350), l = c(60, 40, 20)) + labs(y = "Co-occurrence probability") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_blank(), axis.text.y = element_text(size = 12),
        legend.position = "none", panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  panel_border(remove = F, colour = "black")

ggplot(out, aes(x = pnorm(Z.Score), col = diet.pair)) + geom_density()

# pair diet side-by-side
ggplot(out[out$type != "ZERO",], aes(y = pnorm(Z.Score), x = diet.match, fill = diet.match)) + 
  geom_boxplot(notch = T) + facet_grid(type~time, scales = 'free') + 
  scale_fill_hue(h = c(120, 250, 350), l = c(60, 40, 20)) + labs(y = "Co-occurrence probability", x = "Pair diet") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = "none", panel.background = element_blank(), panel.border = element_rect(fill = NA)) 

out %>% group_by(time) %>% summarise(agg = percpos(Z.Score))
out %>% group_by(diet.match) %>% summarise(agg = percpos(Z.Score))

##### Permutation test code template. ####
outp <- out %>% filter(type == "Segregation", diet.match == TRUE)

d1 <- outp$Z.Score
d1 <- pnorm(d1)
d2 <- outp$time

anv <- permanova(d1, d2, n = 10000)
a <- anova(lm(d1~d2))[["F value"]][1]
length(which(anv> a))/n

permanova <- function(d1, d2, n = 10000) {
  perm <- map(1:n, ~sample(d1, size = length(d1)))
  t <- map(perm, ~lm(resp~indep, data = data.frame(resp = ., indep = d2)))
  anv <- map(t, anova) %>% lapply(`[[`, "F value") %>% sapply(`[`, 1) ### get out F values.
  quantile(anv, 0.99)
  return(anv)
}



#####