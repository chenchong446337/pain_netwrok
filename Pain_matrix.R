## The script is for pain brain area network analysis
## data are from Allen Brian atlas
## created 05282020
## last updated: 05282020 (new way to retrieve data)


library(httr)
library(jsonlite)
library(plyr)
library(tidyverse)
library(igraph)
library(statnet)
library(svMisc)
library(crayon)
library(scales)
library(cowplot)
library(corrplot)



## brain areas which are involved in pain perception
vertex_sym_CH <- c("ACAd","ACAv","SSp", "SSs","MOp", "MOs", "AI", "RSPd","RSPv", # belong to isocortex
                "CA1","CA3","DG", # belong to HPF
                "CLA", "BLA", # belong to cortical subplate
                "CP", "ACB", "AAA", "CEA") ## belong to cerebral nuclei

vertex_sym_BS <- c("VAL", "VM", "VPL", "VPM", "LP", "PO", "AV", "AM", "LD", "MD", "PVT","PT", "MH", "LH", "ZI", "PB","PG", "RN", "RPO", "PAG")               

vertex_sym_CB <- c("CENT", "CUL", "PYR", "NOD", "SIM","AN","PRM", "PFL", "FN", "IP", "IO")

vertex_sym <- c(vertex_sym_CH, vertex_sym_BS, vertex_sym_CB)
vertex_str <- c(rep("CH", length(vertex_sym_CH)), rep("BS", length(vertex_sym_BS)), rep("CB", length(vertex_sym_CB)) )
rm(vertex_sym_CH, vertex_sym_BS, vertex_sym_CB)

## download related data to local drive----
for (i in seq_along(vertex_sym)){
  progress(i, progress.bar = TRUE)
  Sys.sleep (0.01)
  if (i ==length(vertex_sym)) cat("Done!\n")
  
  dat_exp  <- str_c ("http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure[injection_structures$eq", vertex_sym[i],"][transgenic_lines$eq0][primary_structure_only$eqtrue]") %>% 
    GET()
  ID_exp <- fromJSON(rawToChar(dat_exp$content)) %>% 
    .$msg %>% 
    .$id
  
  if (length(ID_exp)==0) {
    cat(str_c("No experiment found for",red(vertex_sym[i]), "in WT mice!\n",  sep = " "))
    dat_exp  <- str_c ("http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure[injection_structures$eq", vertex_sym[i],"][primary_structure_only$eqtrue]") %>% 
      GET()
    ID_exp <- fromJSON(rawToChar(dat_exp$content)) %>% 
      .$msg %>% 
      .$id 
  } 
  if (length(ID_exp)==0) {
    cat(str_c("No experiment found for",red(vertex_sym[i]), "!\n", sep = " "))
    next
    
  } else {
    for (j in seq_along(ID_exp)) {
      dat_sta <- str_c("http://api.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json?criteria=[section_data_set_id$eq", ID_exp[j], "],[is_injection$eqfalse]&num_rows=5000&include=structure") %>% 
        GET()
      value_projection <- fromJSON(rawToChar(dat_sta$content), flatten = T) %>% 
        .$msg %>% 
        as_tibble() 
      write.csv(value_projection, str_c(vertex_sym[i],"_",ID_exp[j],".csv" ))
    }
  }
}


## function to extract the projection values on line----

dat_project_value <- NULL
vertex_coordinate <- vector(mode = "list", length = length(vertex_sym))
for (i in seq_along(vertex_sym)) {
  progress(i, progress.bar = TRUE)
  Sys.sleep (0.01)
  if (i ==length(vertex_sym)) cat("Done!\n")
  
  sym_check <- vertex_sym[-i]
  dat_exp  <- str_c ("http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure[injection_structures$eq", vertex_sym[i],"][transgenic_lines$eq0][primary_structure_only$eqtrue]") %>% 
    GET()
  dat_exp_detail <- fromJSON(rawToChar(dat_exp$content)) %>% 
    .$msg 
  ID_exp <- dat_exp_detail$id
  
  if (length(ID_exp)==0) {
    cat(str_c("No experiment found for",red(vertex_sym[i]), "in WT mice!\n",  sep = " "))
    dat_exp  <- str_c ("http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure[injection_structures$eq", vertex_sym[i],"][primary_structure_only$eqtrue]") %>% 
      GET()
    dat_exp_detail <- fromJSON(rawToChar(dat_exp$content)) %>% 
      .$msg
    ID_exp <- dat_exp_detail$id
    } 
  
  if (length(ID_exp)==0) {
    cat(str_c("No experiment found for",red(vertex_sym[i]), "!\n", sep = " "))
    vertex_coordinate[[i]] <- NA
    next
    
  } else {
    vertex_coordinate[[i]] <- dat_exp_detail %>% 
      filter(id==ID_exp[1]) %>% 
      .$`injection-coordinates` %>% 
      unlist()
   # for (j in seq_along(ID_exp)) {
   #   dat_sta <- str_c("http://api.brain-map.org/api/v2/data/ProjectionStructureUnionize/query.json?criteria=[section_data_set_id$eq", ID_exp[j], "],[is_injection$eqfalse]&num_rows=5000&include=structure") %>% 
   #     GET()
   #   value_projection <- fromJSON(rawToChar(dat_sta$content), flatten = T) %>% 
   #     .$msg %>% 
   #     as_tibble() %>% 
   #     filter(structure.acronym %in% sym_check & hemisphere_id == "2") %>% 
   #     select(structure.acronym, normalized_projection_volume, projection_volume) %>% 
   #     add_column(structure.acronym.IN = vertex_sym[i],.before = "structure.acronym" ) %>% 
   #     add_column(Group = vertex_str[i])
   #   
   #   dat_project_value <- rbind(dat_project_value, value_projection)
   # }
  }
}

# write.csv(dat_project_value, "dat_project_value1.csv")




### start from here with saved data--------

dat_project_value <- NULL
c_select <- function(data){
  data %>% 
    filter(structure.acronym %in% sym_check & hemisphere_id == "2") %>% 
    select(structure.acronym, normalized_projection_volume, projection_volume) 
}
for (i in seq_along(vertex_sym)) {
  progress(i, progress.bar = TRUE)
  Sys.sleep (0.01)
  if (i ==length(vertex_sym)) cat("Done!\n")
  sym_check <- vertex_sym
  
  
  dat_projection <- list.files(pattern = str_c("^",vertex_sym[i], ".*.csv")) %>% 
    sapply(., read.csv, simplify = F) %>% 
    mapply(c_select, ., SIMPLIFY = F) %>% 
    do.call(rbind, .) %>%
    as_tibble() %>% 
    add_column(structure.acronym.IN = vertex_sym[i],.before = "structure.acronym" ) %>% 
    add_column(Group = vertex_str[i]) 
  
  dat_project_value <- rbind(dat_project_value, dat_projection)
  
}

## log10 the value for visualization (Ipilateral)----

dat_project_value_filter<- dat_project_value %>% 
  rename(From = structure.acronym.IN, to = structure.acronym) %>% 
  select(From, to,normalized_projection_volume ) %>% 
  mutate(value = log10(normalized_projection_volume)) %>% 
  mutate(value=replace(value, value < -3.5, -3.5)) %>% 
  mutate(value=replace(value, value > -0.5, -0.5)) %>% 
  mutate(From = factor(From, levels = rev(vertex_sym))) %>% 
  mutate(to=factor(to, levels = vertex_sym)) 

## plot the distribution of projection strength
dat_project_value_hist<- dat_project_value %>% 
  rename(From = structure.acronym.IN, to = structure.acronym) %>% 
  select(From, to,normalized_projection_volume ) %>% 
  filter(From!=to) %>% 
  filter(normalized_projection_volume!=0 & normalized_projection_volume>0.0001 ) %>% 
  ddply(., .(From, to), summarise, mean_value = mean(normalized_projection_volume)) %>% 
  ggplot(., aes(mean_value)) +
  geom_histogram(bins = 100)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


p_plot <- ggplot(dat_project_value_filter, aes(to,From))+
  geom_tile(aes(fill=value)) +
  labs(x="", y = "") +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = c("steelblue4", "yellow", "tomato3"))+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))


## for contralateral projection-----
dat_project_value_contral <- NULL
c_select_contral <- function(data){
  data %>% 
    filter(structure.acronym %in% sym_check & hemisphere_id == "1") %>% 
    select(structure.acronym, normalized_projection_volume, projection_volume) 
}
for (i in seq_along(vertex_sym)) {
  progress(i, progress.bar = TRUE)
  Sys.sleep (0.01)
  if (i ==length(vertex_sym)) cat("Done!\n")
  sym_check <- vertex_sym
  
  
  dat_projection <- list.files(pattern = str_c("^",vertex_sym[i], ".*.csv")) %>% 
    sapply(., read.csv, simplify = F) %>% 
    mapply(c_select_contral, ., SIMPLIFY = F) %>% 
    do.call(rbind, .) %>%
    as_tibble() %>% 
    add_column(structure.acronym.IN = vertex_sym[i],.before = "structure.acronym" ) %>% 
    add_column(Group = vertex_str[i]) 
  
  ## one ID to each experiment
  num_ID <- list.files(pattern = str_c("^",vertex_sym[i], ".*.csv")) %>% 
    length()
  dat_project_value_contral <- rbind(dat_project_value_contral, dat_projection)
  
}
## log10 the value for visualization (contral)

dat_project_value_filter_contral<- dat_project_value_contral %>% 
  rename(From = structure.acronym.IN, to = structure.acronym) %>% 
  select(From, to,normalized_projection_volume ) %>% 
  mutate(value = log10(normalized_projection_volume)) %>% 
  mutate(value=replace(value, value < -3.5, -3.5)) %>% 
  mutate(value=replace(value, value > -0.5, -0.5)) %>% 
  mutate(From = factor(From, levels = rev(vertex_sym))) %>% 
  mutate(to=factor(to, levels = vertex_sym)) 



p_plot_contral <- ggplot(dat_project_value_filter_contral, aes(to,From))+
  geom_tile(aes(fill=value)) +
  scale_x_discrete(position = "top") +
  labs(x="", y = "") +
  scale_fill_gradientn(colours = c("steelblue4", "yellow", "tomato3"))+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

p_heat <- plot_grid(p_plot, p_plot_contral, ncol = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_heat.pdf", width = 300/25.6, height = 120/25.6, family = "Arial")
p_heat
dev.off()

## comparing the projection strength between ipsilateral and contralateral
ratio_strong_connection <- dat_project_value %>% 
  filter(projection_volume> 0.0001) %>% 
  nrow()

ratio_strong_connection_contra <- dat_project_value_contral %>% 
  filter(projection_volume> 0.0001) %>% 
  nrow()

ratio_strong_connection/(ratio_strong_connection + ratio_strong_connection_contra)
## correlation between distance and projection strength-------

## save the coodinate of each vertex
# dat_vertex_coord <- do.call(rbind, vertex_coordinate) %>%
#   dist(., method = "euclidean") %>%
#   cmdscale(., k=2) %>%
#   as_tibble() %>%
#   rename(x=V1, y = V2) %>%
#   add_column(vertex = vertex_sym, .before = "x")
# # 
# write.csv(dat_vertex_coord, "dat_vertex_coord.csv", row.names=FALSE)

dat_vertex_coord <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/R_script/pain_network/dat_vertex_coord.csv")

vertex_dist <- rep(0, nrow(dat_project_value_filter))
for (i in c(1:nrow(dat_project_value_filter))){
  ID_x <- dat_project_value_filter$From[i]
  ID_y <- dat_project_value_filter$to[i]
  if (ID_x ==ID_y){
    vertex_dist[i] = 0 
  } else{
    position_xy <- dat_vertex_coord %>% 
      filter(vertex==ID_x|vertex==ID_y)
    vertex_dist[i] <- as.numeric(dist(position_xy[,c(2,3)], method = "euclidean"))
    
  }
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

dat_project_value_filter1 <- dat_project_value_filter %>% 
  add_column(dis = vertex_dist) %>% 
  add_column(projection_volume = dat_project_value$projection_volume) %>% 
  mutate(projection_volume = ifelse( From == to, 0, projection_volume)) %>% 
  mutate(From=factor(From, levels = vertex_sym)) %>% 
  group_by(From) %>%
  group_modify(~ mutate(., dist_nor = range01(dis)))
  
dat_project_value_filter1_contral <- dat_project_value_filter_contral %>% 
  add_column(dis = vertex_dist) %>% 
  add_column(projection_volume = dat_project_value$projection_volume) %>% 
  mutate(projection_volume = ifelse( From == to, 0, projection_volume)) %>% 
  mutate(From=factor(From, levels = vertex_sym)) %>% 
  group_by(From) %>%
  group_modify(~ mutate(., dist_nor = range01(dis))) %>% 
  mutate(From = factor(From, levels = rev(vertex_sym))) %>% 
  mutate(to=factor(to, levels = vertex_sym)) 
## heat plot of the distance and connection
p_dis_con <- ggplot(dat_project_value_filter1)+
  geom_tile(aes(to,From, fill=dist_nor)) +
  labs(x="", y = "") +
  scale_fill_gradientn(colours = c( "tomato3", "yellow","steelblue4"))+
  geom_point(aes(From, to,size=ifelse(projection_volume >= 0.01, "dot", "no_dot")), shape=1) +
  scale_size_manual(values=c(dot=1, no_dot=NA), guide="none") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_dis_con.pdf", width = 90/25.6, height =90/25.6, family = "Arial")
p_dis_con
dev.off()
## get the correlation of each

vertex_cor <- dat_project_value_filter1 %>% 
  filter(normalized_projection_volume>0.001) %>% 
  filter(From!=to) %>% 
  mutate(y = log10(normalized_projection_volume)) %>% 
  group_map(~ cor(.x$dist_nor, .x$y)) %>% 
  #group_map(~ unlist(cor.test(.x$dist_nor, .x$y )[c("estimate", 'p.value')])) %>% 
  do.call(rbind, .) 
  # as_tibble() %>% 
  # filter(p.value < 0.05)

vertex_cor_contral <- dat_project_value_filter1_contral %>% 
  filter(normalized_projection_volume>0.001) %>% 
  filter(From!=to) %>% 
  mutate(y = log10(normalized_projection_volume)) %>% 
  group_map(~ cor(.x$dist_nor, .x$y)) %>% 
  #group_map(~ unlist(cor.test(.x$dist_nor, .x$y )[c("estimate", 'p.value')])) %>% 
  do.call(rbind, .) 
  #as.tibble() %>% 
  #filter(p.value < 0.05)

dat_cor <- tibble(Group= rep(c("Ipsi", "Contral"), each = length(vertex_cor)), cor=c(vertex_cor, vertex_cor_contral)) %>% 
  mutate(Group = factor(Group, levels = c("Ipsi", "Contral")))

dat_cor_sta <- ddply(dat_cor, .(Group), summarise, n=length(cor),mean=mean(cor),sd=sd(cor),se=sd(cor)/sqrt(length(cor)))
p_cor_dis <- ggplot(dat_cor, aes(Group, cor, colour=Group))+
  geom_violin()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="correlation coefficient (r)")+
  theme_classic()+
  theme(legend.position = 'none')

p_cor_exp <- dat_project_value_filter1 %>% 
  filter(From=="ACAd") %>% 
  filter(normalized_projection_volume>0.0001) %>% 
  ggplot(., aes(dist_nor, log10(normalized_projection_volume)))+
  geom_point()+
  geom_smooth(method = "lm", colour="red")+
  labs(x="Normalized distance between areas", y = "Normalized projection volume (log10)")+
  theme_classic()

p_cor <- plot_grid(p_cor_exp, p_cor_dis, rel_widths = c(1, 0.5))

cor_acad <- dat_project_value_filter1 %>% 
  filter(From=="ACAd") %>% 
  filter(normalized_projection_volume>0.0001) %>% 
  cor.test(~ dist_nor + log10(normalized_projection_volume), data = .)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_cor.pdf", width = 140/25.6, height = 80/25.6, family = "Arial")
p_cor
dev.off()
## correlation analysis of the similarity of target and source input-----
dat_project_value_filter_wide <- dat_project_value_filter %>% 
  select(-normalized_projection_volume) %>% 
  ddply(., .(From, to), summarise, value=mean(value)) %>% 
  pivot_wider(names_from = to, values_from = value) %>% 
  select(vertex_sym)

res <- cor(dat_project_value_filter_wide)
col1<- col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                     "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                     "#4393C3", "#2166AC", "#053061"))(200))
corrplot(res, method = "square",type = "full", tl.col = "black",col=col1 )

dat_project_value_filter_wide1 <- dat_project_value_filter %>% 
  select(-normalized_projection_volume) %>% 
  ddply(., .(From, to), summarise, value=mean(value)) %>% 
  pivot_wider(names_from = From, values_from = value) %>% 
  select(vertex_sym)

res1 <- cor(dat_project_value_filter_wide1)

corrplot(res1, method = "square",type = "full", tl.col = "black", col=col1)

## convert data to network object----
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
dat_project_value1 <- dat_project_value %>% 
  filter(structure.acronym.IN!=structure.acronym) %>% 
  filter(projection_volume > 0.01) %>%  
  mutate(Group = factor(Group, levels = c("CH", "BS", "CB"))) %>% 
  rename(From = structure.acronym.IN, to = structure.acronym) %>% 
  mutate(From = factor(From, levels = vertex_sym)) %>% 
  select(From, to,normalized_projection_volume, Group ) %>% 
  ddply(., .(From, to, Group), summarise, value=mean(normalized_projection_volume)) %>% 
  mutate(color=ifelse(Group=="CH", "brown3", ifelse(Group=="BS", "aquamarine4", "steelblue4")))

  


pain_net <- dat_project_value1[, 1:2] %>% 
  graph.data.frame(., vertices=dat_vertex_coord, directed=TRUE)

## set edge  colors by input and output direction
# x <- t(apply(dat_project_value1[, 1:2], 1, sort))
# x <- paste0(x[, 1], x[, 2])
# y <- x[duplicated(x)]
# COLS <- ifelse(x %in% y, "indianred4", "deepskyblue4")

## set edge color same as the vertex color

E(pain_net)$color <- dat_project_value1$color

E(pain_net)$weight <- dat_project_value1$value

V(pain_net)$Group <- vertex_str[order(match(vertex_sym, names(V(pain_net))))]
vertex_color <- rep(0, length(vertex_str))
vertex_color[which(V(pain_net)$Group=="CH")] ="brown3"
vertex_color[which(V(pain_net)$Group=="BS")] ="aquamarine4"
vertex_color[which(V(pain_net)$Group=="CB")] ="steelblue4"
V(pain_net)$color <- vertex_color

lay_fr <- layout.norm(as.matrix(dat_vertex_coord[,2:3]))

par(mar=c(0,0,0,0)+0.2)
plot(pain_net, layout=lay_fr,vertex.col=V(pain_net)$color, vertex.size=5, edge.width = log10(E(pain_net)$weight)*5, edge.color = adjustcolor(E(pain_net)$color,0.8), edge.arrow.size= 0.2, edge.curved = T )
legend("bottomleft", legend = c("CH", "BS","CB"), col=c("brown3", "aquamarine4", "steelblue4"), pch=19, pt.cex=1.5, bty="n")

par(mar=c(5,4,4,2)+0.1)


# save the figure

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_network.pdf", width = 120/25.6, height = 90/25.6, family = "Arial")

dev.off()

## analyze the correlation between distance and projection

## describe analysis, table of centrality-----
degree_pain_net <- igraph::degree(pain_net, mode = "total")
closness_pain_net <- igraph::closeness(pain_net, mode = "all")
betwness_pain_net <- igraph::betweenness(pain_net, directed = T)

dat_pain_centra <- tibble(Vertex = names(degree_pain_net) ,Degree = degree_pain_net, Closness = closness_pain_net, Betweeness= betwness_pain_net)
## create models and compare with the pain net 

pain_net_rnd <- igraph::erdos.renyi.game(length(V(pain_net)), graph.density(pain_net), type = "gnp",directed = T )

# pain_net_smw <- mean(igraph::degree(pain_net, mode = "out")) %>% 
  igraph::sample_smallworld(1, length(V(pain_net)), nei =., p=0.3)

## Generate small world random network
n <- gorder(pain_net)
m <- gsize(pain_net)

# Work out how theye divide into each other
rem <- m %% n
div <- m %/% n

set.seed(123)
if(rem != 0) {
  pain_net_smw <- sample_smallworld(1, n, div+1, p = 0.05)
  # Randomly delete the unwanted edges. Should be quite homegenous
  pain_net_smw <- delete_edges(pain_net_smw, sample(1:gsize(pain_net_smw), size = gsize(pain_net_smw) - m))
} else {
  pain_net_smw <- sample_smallworld(1, n, div, p = 0.001)
}

gorder(pain_net_smw)
gsize(pain_net_smw)

par(mar=c(0,0,0,0)+0.2)
plot(pain_net_smw, layout=lay_fr,vertex.col=V(pain_net)$color, vertex.size=5, edge.width = 0.2, edge.color = adjustcolor(E(pain_net)$color,0.8), edge.arrow.size= 0.2, edge.curved = T )

par(mar=c(5,4,4,2)+0.1)
## generate scale free graph with same nodes and edges
genOutSeq <- function(n, m) {
  n <- n-1 # Shift it along
  rem <- m %% n
  c(0, rep(m%/%n + 1, rem), rep(m%/%n, n - rem))
}

set.seed(11)
pain_net_sfg <- sample_pa(gorder(pain_net), power = 0.5,
                          out.seq = genOutSeq(gorder(pain_net), gsize(pain_net)),
                          algorithm = "psumtree-multiple", directed = T)

V(pain_net_sfg)$name <- V(pain_net)$name

gorder(pain_net_sfg)
gsize(pain_net_sfg)
par(mar=c(0,0,0,0)+0.2)
plot(pain_net_sfg, vertex.col=V(pain_net)$color, vertex.size=5, edge.width = 0.2,  edge.arrow.size= 0.2, edge.curved = T )

par(mar=c(5,4,4,2)+0.1)

cluster_walktrap(pain_net_sfg)
mean_distance(pain_net_sfg, directed = T, unconnected = F)
hist(igraph::degree(pain_net_sfg))


# pain_net_sfg <- igraph::sample_pa(length(V(pain_net)), out.dist = c(0.25, 0.5, 0.25),directed = T )

## calculate the reciprocity of graphs
recip_pain_net <- igraph::reciprocity(pain_net)
recip_pain_net_rnd <- igraph::reciprocity(pain_net_rnd)
recip_pain_net_smw <- igraph::reciprocity(pain_net_smw)
recip_pain_net_sfg <- igraph::reciprocity(pain_net_sfg)

## histogram plot of the degree, compare with other three group

d_pain_net <- tibble::enframe(igraph::degree(pain_net, mode = "total")) %>% 
  rename(vertex_sym = name,d_pain = value)
d_pain_net_rnd <- tibble::enframe(igraph::degree(pain_net_rnd, mode = "total")) %>% 
  rename(d_rnd = value)
d_pain_net_smw <- tibble::enframe(igraph::degree(pain_net_smw, mode = "total")) %>% 
  rename(d_smw = value)
d_pain_net_sfg <- tibble::enframe(igraph::degree(pain_net_sfg, mode = "total")) %>% 
  rename(d_sfg = value) 

##cummulative plot
dat_d_pain <- cbind(d_pain_net, d_pain_net_rnd[,2], d_pain_net_smw[,2], d_pain_net_sfg[,2]) %>% 
  pivot_longer(-vertex_sym) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()

dat_d_pain <- cbind(d_pain_net, d_pain_net_rnd[,2], d_pain_net_smw[,2], d_pain_net_sfg[,2]) %>% 
  pivot_longer(-vertex_sym) %>% 
  ggplot(., aes(value,  color=name)) + 
  geom_density()+
  theme_classic()

## compare the concentration of degree----

lorenz = function(d) {
  d = sort(d, dec=TRUE)
  n = length(d);
  s = sum(d);
  p = seq(0, 1, 0.01);
  c = NULL;
  items = NULL;
  i = 1;
  for (x in p) {
    if (x == 0)  {
      items[i] = 0;
      c[i] = 0;
    } else { 
      items[i] = floor(n * x);
      c[i] = sum(d[1:items[i]]) / s;
    }  
    i = i + 1;
  }
  names(c) = p
  list = list(p = p, c = c, items = items); 
}

dat_r <- sapply(list(pain_net, pain_net_rnd, pain_net_sfg, pain_net_smw), igraph::degree, mode = "total") %>% 
  apply(., 2, lorenz)

p_dat_r <- cbind(dat_r[[1]]$p, dat_r[[1]]$c, dat_r[[2]]$c,dat_r[[3]]$c,dat_r[[4]]$c) %>% 
  as_tibble() %>% 
  rename(Per = V1, Pain_net = V2, Pain_rnd = V3, Pain_sfg = V4, Pain_smw = V5) %>% 
  pivot_longer(!Per, names_to = 'Group', values_to = "Value") %>% 
  ggplot(., aes(Per, Value, colour = Group))+
  geom_line()+
  theme_classic()+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_dat_r.pdf", width = 75/25.6, height = 75/25.6, family = "Arial")
p_dat_r
dev.off()

## compare different models----
# for degree
degree_compare <- sapply(list(pain_net, pain_net_rnd, pain_net_sfg, pain_net_smw), igraph::degree, mode = "total") %>% 
                       apply(., 2, mean)    

transitivity_compare <- sapply(list(pain_net, pain_net_rnd, pain_net_sfg, pain_net_smw), igraph::transitivity, type = "local", isolates="NaN") %>% 
  apply(., 2, mean, na.rm = T) 

density_compare <- sapply(list(pain_net, pain_net_rnd, pain_net_sfg, pain_net_smw), igraph::edge_density, loops = FALSE)

reciprocity_compare <- sapply(list(pain_net, pain_net_rnd, pain_net_sfg, pain_net_smw), igraph::reciprocity)

## histogram plot of the clustering coefficient, compare with other three group

c_pain_net <- tibble::enframe(igraph::transitivity(pain_net,type = "local", isolates="NaN")) %>% 
  rename(vertex_sym = name,c_pain = value)
c_pain_net_rnd <- tibble::enframe(igraph::transitivity(pain_net_rnd,type = "local", isolates="NaN")) %>% 
  rename(c_rnd = value)
c_pain_net_smw <- tibble::enframe(igraph::transitivity(pain_net_smw,type = "local", isolates="NaN")) %>% 
  rename(c_smw = value)
c_pain_net_sfg <- tibble::enframe(igraph::transitivity(pain_net_sfg,type = "local", isolates="NaN")) %>% 
  rename(c_sfg = value) 

##cummulative plot

dat_c_pain <- cbind(c_pain_net, c_pain_net_rnd[,2], c_pain_net_smw[,2], c_pain_net_sfg[,2]) %>% 
  pivot_longer(-vertex_sym) %>% 
  filter(value!='NaN' & value!=0) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()

dat_c_pain <- cbind(c_pain_net, c_pain_net_rnd[,2], c_pain_net_smw[,2], c_pain_net_sfg[,2]) %>% 
  pivot_longer(-vertex_sym) %>% 
  filter(value!='NaN' & value!=0) %>% 
  ggplot(., aes(value, color=name)) + 
  geom_density()+
  theme_classic()

p_hist_net <- plot_grid(dat_d_pain, dat_c_pain, ncol = 2)

# save the figure
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("p_hist_net.pdf", width = 210/25.6, height = 80/25.6, family = "Arial")
p_hist_net
dev.off()

## find the hubs----
dat_hub <- hub_score(pain_net, scale = TRUE, weights = NULL,options = arpack_defaults)$vector %>% 
  enframe(., name = "name", value = "value") %>% 
  arrange(desc(value))

dat_hub_sfg <- hub_score(pain_net_sfg, scale = TRUE, weights = NULL,options = arpack_defaults)$vector %>% 
  enframe(., name = "name", value = "value") %>% 
  arrange(desc(value))
## interactive plot
tkplot(pain_net, layout=lay_fr,vertex.col=V(pain_net)$color, vertex.size=5, edge.width = log10(E(pain_net)$weight)*5, edge.color = adjustcolor(E(pain_net)$color,0.8), edge.arrow.size= 0.2, edge.curved = T )



## network description
edge_density(pain_net, loops = F)
igraph::diameter(pain_net, directed = T, weights = NA)
transitivity(pain_net, type = "average")
transitivity(pain_net_rnd, type = "average")

edge_density(pain_net_rnd, loops = F)
## the proportion of reciprocated ties
reciprocity(pain_net)

dyad_census(pain_net)
dyad_census(pain_net_rnd)
triad_census(pain_net)
V(pain_net)$size<-igraph::degree(pain_net, mode="all", normalized = T)
d_pain_net_in <- tibble::enframe(igraph::degree(pain_net, mode = "in")) %>% 
  rename(vertex_sym = name,d_in = value)
d_pain_net_out <- tibble::enframe(igraph::degree(pain_net, mode = "out")) %>% 
  rename(d_out = value)
d_pain_net <- tibble::enframe(igraph::degree(pain_net, mode = "all")) %>% 
  rename(d_total = value)


dat_d_pain <- cbind(d_pain_net_in, d_pain_net_out[,2], d_pain_net[,2]) %>% 
  pivot_longer(-vertex_sym) %>% 
  ggplot(., aes(value, color=name)) + stat_ecdf(geom = "step")

## random model
pain_net_ran <- erdos.renyi.game(length(V(pain_net)), length(E(pain_net)), type="gnm", directed = T)

plot(d_pain_net, type = "l")
igraph::closeness(pain_net)
igraph::betweenness(pain_net)
igraph::evcent(pain_net)
igraph::centralization.betweenness(pain_net)
igraph::centralization.degree(pain_net)
igraph::components(pain_net, mode = "weak")
igraph::transitivity(pain_net, type = "local", isolates="zero")
sum(igraph::degree(pain_net, mode="all"))

cluster_edge_betweenness(pain_net, weights =NULL,
                         directed = TRUE, edge.betweenness = TRUE, merges = TRUE,
                         bridges = TRUE, modularity = TRUE, membership = TRUE)
igraph::clique.number(pain_net)
igraph::cliques(pain_net, min=3)

neighbors(pain_net, "PG", mode = "in")
## community detection---
ceb <- cluster_edge_betweenness(pain_net, directed = T) 
dendPlot(ceb, mode="hclust")

clp <- cluster_label_prop(pain_net) 
plot(clp, pain_net, )

cw_pain <- cluster_walktrap(pain_net)
modularity(cw_pain)
vertex_memb <- membership(cw_pain) %>% 
  as.vector() 

dat_vertex_member <- tibble(Vertex = vertex_sym, Comminity = vertex_memb)

par(mar=c(0,0,0,0)+0.2)
plot(pain_net1, layout=lay_fr,vertex.col=unname(membership(cw_pain)), vertex.size=10, edge.arrow.size= 0.2,  edge.curved = T)
legend("bottomleft", legend = c("CH", "BS","CB"), col=c("brown3", "aquamarine4", "steelblue4"), pch=19, pt.cex=1.5, bty="n")

par(mar=c(5,4,4,2)+0.1)

## k-cores analysis---
k_pain_net <- subgraph.edges(pain_net, E(pain_net)[weight>0.3])
pain_coreness <- graph.coreness(k_pain_net)
V(k_pain_net)$color <- pain_coreness +1

plot(k_pain_net, layout=lay_fr,vertex.col=membership(cw_pain), vertex.size=10, edge.arrow.size= 0.2,  edge.curved = T)

## 
pain_net_ran <- erdos.renyi.game(length(V(pain_net)), length(E(pain_net)), type = "gnm", directed = T)
pain_net_ran <- erdos.renyi.game(11, 0.81, type = "gnp", directed = T)

edge_density(pain_net_ran, loops = F)

igraph::graph.density(pain_net_ran)
mean(igraph::degree(pain_net))
mean(igraph::degree(pain_net_ran))
igraph::transitivity(pain_net_ran)

## cutpoint of the network
articulation_points(pain_net)
clique.number(pain_net)

## calculate the assortivity----
assortativity_degree(pain_net, directed = TRUE)
assortativity_degree(pain_net_rnd, directed = TRUE)
assortativity_degree(pain_net_sfg, directed = TRUE)
assortativity_degree(pain_net_smw, directed = TRUE)

assortativity_nominal(pain_net,as.numeric(as.factor(V(pain_net)$Group)),directed = TRUE)

## remove some vertices
pain_net_rm <- delete_vertices(pain_net, c("CENT", "CUL","PYR", "NOD", "SIM", "AN", "PRM", "PFL"))
assortativity_degree(pain_net_rm, directed = TRUE)

library(Matrix)
A <- get.adjacency(pain_net)
image(A)

i <- order( igraph::degree(pain_net) )
image(A[i,i])

## test how remove node change the diameter----

n_rm <- ceiling(seq(from = 0, to = 0.2, 0.02)*length(V(pain_net)))

dis_mean_total1 <- rep(0, 100*length(n_rm))
size_ratio_total1 <- rep(0, 100*length(n_rm))
aver_ratio_total1 <- rep(0, 100*length(n_rm))

for (i in c(1: length(n_rm))){
  for (j in c(1:100)){
    vertex_rm <- sample(as_ids(V(pain_net)), n_rm[i])
    p_net <- delete_vertices(pain_net, vertex_rm) 
    dis_mean <- mean_distance(p_net, directed = T, unconnected = F)
    dis_mean_total1[(i-1)*100 + j] <- dis_mean
    
    ## copare the cluster size change
    max_pain_net <- cluster_walktrap(pain_net) %>% 
      sizes() %>% 
      max()/gorder(pain_net)
    
    size_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      max()/gorder(p_net)/max_pain_net
    size_ratio_total1[(i-1)*100 + j] <- size_ratio
    
    ## average of isolated nodes
    aver_pain_net <-cluster_walktrap(pain_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()
    
    average_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()/aver_pain_net
    
    aver_ratio_total1[(i-1)*100 + j] <- average_ratio
  }
}

## rm by hubs
rm_hub <- igraph::degree(pain_net) %>% 
  .[order(., decreasing = TRUE)] %>% 
  names()

n_rm <- ceiling(seq(from = 0, to = 0.2, 0.02)*length(V(pain_net)))

dis_mean_total2 <- rep(0, 100*length(n_rm))
size_ratio_total2 <- rep(0, 100*length(n_rm))
aver_ratio_total2 <- rep(0, 100*length(n_rm))

for (i in seq_along(n_rm)){
  for (j in c(1:100)){
    if (n_rm[i] == 0){
      p_net <- pain_net
    } else{
      vertex_rm <- rm_hub[1:i]
      p_net<- delete_vertices(pain_net, vertex_rm) 
    }
    
    dis_mean <- mean_distance(p_net, directed = T, unconnected = F)
    dis_mean_total2[(i-1)*100 + j] <- dis_mean
    
    ## copare the cluster size change
    max_pain_net <- cluster_walktrap(pain_net) %>% 
      sizes() %>% 
      max()/gorder(pain_net)
    
    size_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      max()/gorder(p_net)/max_pain_net
    size_ratio_total2[(i-1)*100 + j] <- size_ratio
    
    ## average of isolated nodes
    aver_pain_net <-cluster_walktrap(pain_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()
    
    average_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()/aver_pain_net
      
    aver_ratio_total2[(i-1)*100 + j] <- average_ratio
  }
}


f_group <- rep(rep(seq(from = 0, to = 0.2, 0.02), each = 100), 2)
dis_mean_total <- c(dis_mean_total1, dis_mean_total2)
size_ratio_total <- c(size_ratio_total1, size_ratio_total2)
aver_ratio_total <- c(aver_ratio_total1, aver_ratio_total2)

dat_mean_dis <- tibble(Group = f_group, dis = dis_mean_total) %>% 
  mutate(Group = as.factor(Group)) %>% 
  mutate(Group1 = rep(c("Random", "Attack"), each= 1100)) %>% 
  ddply(., .(Group, Group1), summarise,n=length(dis),mean=mean(dis),sd=sd(dis),se=sd(dis)/sqrt(length(dis))) %>% 
  ggplot(., aes(Group, mean, group = Group1, colour = Group1))+
  geom_line(color="grey")+
  geom_point()+
  labs(x="f", y="Mean distance")+
  theme_classic()+
  scale_y_continuous(limits = c(7, 12))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("dat_mean_dis.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
dat_mean_dis
dev.off()
## plot the size ratio
dat_mean_fra <- tibble(Group = f_group, f_size = size_ratio_total, aver_size = aver_ratio_total) %>% 
  mutate(Group = as.factor(Group)) %>% 
  mutate(Group1 = rep(c("Random", "Attack"), each= 1100)) %>% 
  pivot_longer(!c(Group, Group1), names_to = "Variable") %>% 
  ddply(., .(Group, Group1, Variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Group, mean, group = interaction(Group1, Variable)))+
  geom_line(color="grey")+
  geom_point(aes(colour = Group1, shape = Variable ))+
  labs(x="f", y="<s> and S")+
  theme_classic()+
  scale_y_continuous(limits = c(0, 2))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("dat_mean_fra.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
dat_mean_fra
dev.off()




## compare with the scale free model
n_rm <- ceiling(seq(from = 0, to = 0.2, 0.02)*length(V(pain_net_sfg)))

dis_mean_total3 <- rep(0, 100*length(n_rm))
size_ratio_total3 <- rep(0, 100*length(n_rm))
aver_ratio_total3 <- rep(0, 100*length(n_rm))

for (i in c(1: length(n_rm))){
  for (j in c(1:100)){
    vertex_rm <- sample(as_ids(V(pain_net_sfg)), n_rm[i])
    p_net <- delete_vertices(pain_net_sfg, vertex_rm) 
    dis_mean <- mean_distance(p_net, directed = T, unconnected = F)
    dis_mean_total3[(i-1)*100 + j] <- dis_mean
    
    ## copare the cluster size change
    max_pain_net <- cluster_walktrap(pain_net_sfg) %>% 
      sizes() %>% 
      max()/gorder(pain_net_sfg)
    
    size_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      max()/gorder(p_net)/max_pain_net
    size_ratio_total3[(i-1)*100 + j] <- size_ratio
    
    ## average of isolated nodes
    aver_pain_net <-cluster_walktrap(pain_net_sfg) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()
    
    average_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()/aver_pain_net
    
    aver_ratio_total3[(i-1)*100 + j] <- average_ratio
  }
}

## rm by hubs
rm_hub <- igraph::degree(pain_net_sfg) %>% 
  .[order(., decreasing = TRUE)] %>% 
  names()

n_rm <- ceiling(seq(from = 0, to = 0.2, 0.02)*length(V(pain_net_sfg)))

dis_mean_total4 <- rep(0, 100*length(n_rm))
size_ratio_total4 <- rep(0, 100*length(n_rm))
aver_ratio_total4 <- rep(0, 100*length(n_rm))

for (i in seq_along(n_rm)){
  for (j in c(1:100)){
    if (n_rm[i] == 0){
      p_net <- pain_net_sfg
    } else{
      vertex_rm <- rm_hub[1:i]
      p_net<- delete_vertices(pain_net_sfg, vertex_rm) 
    }
    
    dis_mean <- mean_distance(p_net, directed = T, unconnected = F)
    dis_mean_total4[(i-1)*100 + j] <- dis_mean
    
    ## copare the cluster size change
    max_pain_net <- cluster_walktrap(pain_net_sfg) %>% 
      sizes() %>% 
      max()/gorder(pain_net_sfg)
    
    size_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      max()/gorder(p_net)/max_pain_net
    size_ratio_total4[(i-1)*100 + j] <- size_ratio
    
    ## average of isolated nodes
    aver_pain_net <-cluster_walktrap(pain_net_sfg) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()
    
    average_ratio <- cluster_walktrap(p_net) %>% 
      sizes() %>% 
      .[. != max(.)] %>% 
      mean()/aver_pain_net
    
    aver_ratio_total4[(i-1)*100 + j] <- average_ratio
  }
}


f_group <- rep(rep(seq(from = 0, to = 0.2, 0.02), each = 100), 4)
dis_mean_total <- c(dis_mean_total1, dis_mean_total3, dis_mean_total2,dis_mean_total4)


dat_mean_dis <- tibble(Group = f_group, dis = dis_mean_total) %>% 
  mutate(Group = as.factor(Group)) %>% 
  mutate(Group1 = rep(c("Random", "Attack"), each= 2200)) %>% 
  mutate(Group2 = c(rep("Pain", 1100), rep("Scale_free", 1100), rep("Pain", 1100), rep("Scale_free", 1100))) %>% 
  ddply(., .(Group, Group1, Group2), summarise,n=length(dis),mean=mean(dis),sd=sd(dis),se=sd(dis)/sqrt(length(dis))) %>% 
  ggplot(., aes(Group, mean, group = interaction(Group1, Group2)))+
  #geom_line(color="grey")+
  geom_point(aes(colour= Group1, shape = Group2))+
  labs(x="f", y="Mean distance")+
  theme_classic()+
  scale_y_continuous(limits = c(7, 12))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))


## for random--
f_group <- rep(rep(seq(from = 0, to = 0.2, 0.02), each = 100), 4)
size_ratio_total <- c(size_ratio_total1,size_ratio_total2, size_ratio_total3, size_ratio_total4)
aver_ratio_total <- c(aver_ratio_total1, aver_ratio_total2, aver_ratio_total3, aver_ratio_total4)

dat_fra_aver <- tibble(Group = f_group, aver_size = aver_ratio_total) %>% 
  mutate(Group = as.factor(Group)) %>% 
  mutate(Group1 = rep(c("Pain", "Scale_free"), each= 1100*2)) %>% 
  mutate(Group2 = c(rep("Random", 1100), rep("Attack", 1100), rep("Random", 1100), rep("Attack", 1100))) %>% 
  ddply(., .(Group, Group1, Group2), summarise,n=length(aver_size),mean=mean(aver_size, na.rm = T),sd=sd(aver_size),se=sd(aver_size)/sqrt(length(aver_size))) %>% 
  ggplot(., aes(Group, mean, group = interaction(Group1, Group2)))+
  geom_line(color="grey")+
  geom_point(aes(color = Group2, shape = Group1))+
  labs(x="f", y="<s>")+
  theme_classic()+
  scale_y_continuous(limits = c(0, 2))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("dat_fra_aver.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
dat_fra_aver
dev.off()

dat_fra_size <- tibble(Group = f_group, f_size = size_ratio_total) %>% 
  mutate(Group = as.factor(Group)) %>% 
  mutate(Group1 = rep(c("Pain", "Scale_free"), each= 1100*2)) %>% 
  mutate(Group2 = c(rep("Random", 1100), rep("Attack", 1100), rep("Random", 1100), rep("Attack", 1100))) %>% 
  ddply(., .(Group, Group1, Group2), summarise,n=length(f_size),mean=mean(f_size, na.rm = T),sd=sd(f_size),se=sd(f_size)/sqrt(length(f_size))) %>% 
  ggplot(., aes(Group, mean, group = interaction(Group1, Group2)))+
  geom_line(color="grey")+
  geom_point(aes(color = Group2, shape = Group1))+
  labs(x="f", y="S")+
  theme_classic()+
  scale_y_continuous()+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/pain_network/")
cairo_pdf("dat_fra_size.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
dat_fra_size
dev.off()

## save data to plot in gephi------
pain_net_edge <- cbind( get.edgelist(pain_net) , round( E(pain_net)$weight, 3 )) %>% 
  as_tibble() %>% 
  rename(Source = V1, Target = V2, Weight = V3)
  
pain_net_node <- tibble(ID = c(1:vcount(pain_net)), NAME = V(pain_net)$name)

setwd("~cchen2/Documents/neuroscience/Pn\ project/R_script/pain_network/gexl/")
write.csv2(pain_net_edge,row.names = F,'pain_net_edge.csv')

write.csv2(pain_net_node,row.names = F,'pain_net_node.csv')


