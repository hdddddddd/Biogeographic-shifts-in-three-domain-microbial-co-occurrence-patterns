library(dplyr)
myTable=read.csv("inter_sparcc0.75_cytoscape.csv")
nodesList=read.csv("Node_list.csv")

newTable = myTable
colnames(newTable) = c("Target","Source","corrWeight","Type")
newTable = cbind(newTable$Source,newTable$Target,
newTable$corrWeight,
as.character(newTable$Type))
colnames(newTable) = c("Source","Target","corrWeight","Type")
newTable = as.data.frame(newTable, stringsAsFactors=FALSE)

myTable2 = rbind(newTable,myTable)
#Combine nodesList to the myTable2 to make it easier subset the table
Source = nodesList %>% as.data.frame()
colnames(Source) = c("Source","SourceCommunity")
myTable2 = merge(myTable2,Source)
Target = nodesList %>% as.data.frame()
colnames(Target) = c("Target","TargetCommunity")
myTable2 = merge(myTable2,Target)
colnames(myTable2) = c("Target","Source","Weight",
"Type","SourceCommunity","TargetCommunity")
myTable = myTable2

# Before Zi-Pi plot calculations, a little variable clean up
rm(myTable2,Source, Target, newTable)
# Create Zi-Pi table
final_table = matrix(NA, ncol = 3, nrow = length(nodesList$Id)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(final_table) = c("Id","Zi","Pi")
myCommunities = unique(myTable$SourceCommunity) %>% as.character(.)
#Zi calculations
for (x in (1:length(myCommunities))){
## Leave only edges coming nodes from the same community
myComm = myCommunities[x] %>% as.character(.)
currentComm = subset(myTable, SourceCommunity==myComm)
currentComm = subset(currentComm, TargetCommunity==myComm)
if (dim(currentComm)[1] == 0){next} #Break out if the current Community has only negative edges to other nodes, or if the node doesn't have any edges to nodes within the same community
uniq_comm_ids = unique(currentComm$Source)
comm_zi = matrix(NA, ncol = 5, nrow = length(uniq_comm_ids)) %>% as.data.frame(.,stringsAsFactors=FALSE)
colnames(comm_zi) = c("Id","NumOfConnections","average","STDev","Zi")
comm_zi$Id = uniq_comm_ids
for (xx in uniq_comm_ids){
myNum = subset(currentComm, Source==xx)
myNum = dim(myNum)
myNum = myNum[1]
comm_zi = within(comm_zi, NumOfConnections[Id == xx] <- myNum)
}
comm_zi_STDEV = sd(comm_zi$NumOfConnections)
if (comm_zi_STDEV == 0){
comm_zi$Zi = NA
} else {
comm_zi$STDev = comm_zi_STDEV
comm_zi$average = ave(comm_zi$NumOfConnections)
comm_zi$Zi = (comm_zi$NumOfConnections - comm_zi$average)/comm_zi$STDev
}
comm_zi$NumOfConnections = NULL
comm_zi$average = NULL
comm_zi$STDev = NULL
final_table = merge(final_table,comm_zi, by="Id", all.x = TRUE, all.y = TRUE)
myValues = coalesce(as.double(final_table$Zi.y),as.double(final_table$Zi.x))
final_table$Zi.x = NULL
colnames(final_table) = c("Id","Pi","Zi")
final_table$Zi = myValues
}
#Clean up after Zi calculations
rm(comm_zi,comm_zi_STDEV,currentComm,myComm,myCommunities,myNum,myValues,uniq_comm_ids,x,xx)
#Pi calculations
for (x in nodesList$Id){
piTotalTable = subset(myTable, Source==x)
total = dim(piTotalTable)[1]
sigma = 0
for (y in unique(piTotalTable$TargetCommunity)){
commConnections = subset(piTotalTable, TargetCommunity==y) %>% dim(.) %>% .[1]
fraction = (commConnections/total)^2
sigma = sigma + fraction
}
featurePiScore = 1-sigma
final_table = within(final_table, Pi[Id == x] <- featurePiScore)
}
inter_zi_pi = merge(final_table,nodesList)

library(ggplot2)
zi_pi <- na.omit(inter_zi_pi)
zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$Zi < 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$Zi > 2.5 & zi_pi$Pi > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(Pi, Zi)) +
geom_point(aes(color = type), alpha = 0.6, size = 2.3) +
scale_color_manual(values = c('gray','#5A9367','#4ea1d3','#e97f02'), limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
theme_bw()+theme(panel.grid = element_blank(), panel.background = element_blank(), legend.key = element_blank())+
labs(x = 'Among-module connectivities (Pi)', y = 'Within-module connectivities (Zi)', color = '') +
geom_vline(xintercept = 0.62,linetype="dotted",color='#6b6b6b') +
geom_hline(yintercept = 2.5,linetype="dotted",color='#6b6b6b')+
ylim(-1,5)+xlim(0,1.3)+
theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10) )
write.csv(zi_pi,"zi_pi.csv")
ggsave('zi_pi.pdf',width = 15,height = 10,units='cm',dpi=300)
