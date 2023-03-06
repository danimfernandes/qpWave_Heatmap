# Batch read all qpWave output files and retrieve rank=0 p-values
rm(list = ls())
file.remove(list.files(path=getwd(),pattern="NA_QPWAVE_OUT$"))
mylist1=list.files(path=getwd(),pattern="_QPWAVE_OUT$")
testList = c()
for(i in list.files(path=getwd(),pattern="_QPWAVE_OUT$")) {
  testList = c(testList,strsplit(i,"_PLUS_")[[1]][1])
  print(strsplit(i,"_PLUS_")[[1]][1])
  print(paste0("    ",strsplit(i,"_PLUS_")[[1]][2]))
}
sampleList = unique(testList)
amodels_table = data.frame(matrix(nrow = 0,ncol = 5))
colnames(amodels_table)=c("A","B","Right","P-value_R0","SNPs")

counter = 1
for(qpwavefile in mylist1) {
  qpwaveload = read.table(file = qpwavefile,fill = T, sep="\n",stringsAsFactors = F)
  lineL = which(grepl(pattern = "left pops:",x = qpwaveload$V1) == TRUE)+1 ##First left population
  lineR = which(grepl(pattern = "right pops:",x = qpwaveload$V1) == TRUE)+1 ##First right individual
  left_pops = qpwaveload[seq(lineL,lineL+1),]
  
  lineRlast = (which(grepl(pattern = " 0 ",x = qpwaveload$V1) == TRUE)-1)[1]
  right_pops = qpwaveload[seq(lineR,lineRlast),]
  
  line0 = which(grepl(pattern = "f4rank: 0",x = qpwaveload$V1) == TRUE)
  pvalue = strsplit(x = as.character(qpwaveload[line0,]),split = " ")[[1]][strsplit(x = as.character(qpwaveload[line0,]),split = " ")[[1]] != ""][8]
  
  linesnp = which(grepl(pattern = "numsnps used:",x = qpwaveload$V1) == TRUE)
  numspns = strsplit(x = as.character(qpwaveload[linesnp,]),split = " ")[[1]][strsplit(x = as.character(qpwaveload[linesnp,]),split = " ")[[1]] != ""][3]
  
  line=paste0(c(left_pops,paste0(right_pops,collapse = ","),pvalue,numspns))
  amodels_table[counter,] = line
  counter = counter + 1
}
  
## Sorting by column index [1] then [2]
amodels_table = amodels_table[order( amodels_table[,1], amodels_table[,2] ),]
amodels_table$`P-value_R0` = as.numeric(amodels_table$`P-value_R0`)
## Replacing values smaller than 1E-299 with 0 to avoid incompatibilities with Excel when manually reordering the heatmap further down
amodels_table[which(amodels_table$`P-value_R0` < 1e-299 & amodels_table$`P-value_R0` > 0),4] = 0
write.table(x = amodels_table,file = "table_batch_qpWave",sep = "\t",quote = F,row.names = F,col.names = T)

# Load "table_batch_qpWave" to create a heatmap
rm(list = ls())

qpWave_batch_load = read.table(file = "table_batch_qpWave", fill= T, sep="",stringsAsFactors = F,header = T)
qpWave_batch_load[,6] = paste0(qpWave_batch_load$A,"---",qpWave_batch_load$B)
qpWave_batch_load[,7] = paste0(qpWave_batch_load$B,"---",qpWave_batch_load$A)

pops1 = unique(qpWave_batch_load$A)
pops2 = unique(qpWave_batch_load$B)

uniquePops = unique(c(pops1,pops2))
uniquePops = sort(uniquePops)

rm1 = which(duplicated(qpWave_batch_load$V6) == TRUE)
rm2 = which(duplicated(qpWave_batch_load$V7) == TRUE)
if(length(c(rm1,rm2)) != 0) {
  qpWave_batch_load = qpWave_batch_load[-c(rm1,rm2),]}

## Heatmap matrix of p-values
qpWave_heat = matrix(nrow=length(uniquePops),ncol=length(uniquePops))
rownames(qpWave_heat) = uniquePops
colnames(qpWave_heat) = uniquePops

sequPops = matrix(ncol=2,nrow=length(uniquePops))
sequPops[,1] = uniquePops
sequPops[,2] = seq(1:length(uniquePops))

combis = combn(uniquePops,2)
it=1
for(colu in seq(1:length(combis[1,]))) {
  if(length(which(qpWave_batch_load$V6 %in% paste0(combis[,colu][1],"---",combis[,colu][2]))) != 0) {
    where = which(qpWave_batch_load$V6 %in% paste0(combis[,colu][1],"---",combis[,colu][2]))
  } else if(length(which(qpWave_batch_load$V7 %in% paste0(combis[,colu][1],"---",combis[,colu][2]))) != 0) {
    where = which(qpWave_batch_load$V7 %in% paste0(combis[,colu][1],"---",combis[,colu][2]))
  }
  pvalue = qpWave_batch_load[where,4]

  popsHere1 = qpWave_batch_load[where,1]
  popsHere2 = qpWave_batch_load[where,2]
  qpWave_heat[popsHere1,popsHere2] = pvalue
  qpWave_heat[popsHere2,popsHere1] = pvalue
  it=it+1
}
## Output the p-value matrix if needed to be ordered or edited outside of R
write.csv(x = qpWave_heat,file = "TESTtable_qpWave_heatmap.csv",quote = F)

## (Un)comment next line to read edited file, allowing to easily import modified labels and order
qpWave_heat = read.csv(file = "table_qpWave_heatmap.csv", stringsAsFactors = F, row.names = 1)
uniquePops = rownames(qpWave_heat)
qpWave_heat = as.matrix(qpWave_heat)
qpWave_heat_col = qpWave_heat
qpWave_heat_col[qpWave_heat_col < 0.01000] = "grey85"
qpWave_heat_col[qpWave_heat_col >= 0.01000 & qpWave_heat_col <0.05] = "#d8d78f"
qpWave_heat_col[qpWave_heat_col >= 0.05000 & qpWave_heat_col <0.10] = "#ceb43d"
qpWave_heat_col[qpWave_heat_col >= 0.1000 & qpWave_heat_col <0.40] = "#c34818"
qpWave_heat_col[qpWave_heat_col >= 0.4000 & qpWave_heat_col <0.70] = "#7b2300"
qpWave_heat_col[qpWave_heat_col >= 0.7000 & qpWave_heat_col <=1] = "#421300"
qpWave_heat_col[is.na(qpWave_heat_col)] = "white"


###  Heatmap matrix with color values
## Change margin sizes here by replacing the 1st and 4th value in "mar"
par(mar = c(13,0.5,0.5,13), mgp = c(3,0.5,0), xpd=NA)

plot(x=0,y=0,type="n",xlab="",ylab="",xlim = c(0,length(qpWave_heat_col[,1])), ylim = c(length(qpWave_heat_col[,1]),0), bty="n", axes = F, xaxs = "i", yaxs = "i")

for(coli in seq(1:length(qpWave_heat_col[1,]))) {
  ## In coli 1, do X rows, etc
  for(i in seq(1:length(qpWave_heat_col[,1]))) {
    rect(xleft = coli-1,ybottom = i-1,xright = coli,ytop = i,col=qpWave_heat_col[i,coli], lwd=0, border = qpWave_heat_col[i,coli]) # if you want borders, change arg border to for example "white"
  }
}

axis(at = seq(1:length(qpWave_heat_col[1,]))-0.5,labels = colnames(qpWave_heat_col),side = 1, tick = F, las = 2, cex.axis=0.8)
axis(at = seq(1:length(qpWave_heat_col[1,]))-0.5,labels = rownames(qpWave_heat_col),side = 4, tick = F, las = 2, cex.axis=0.8)

## Add legend to bottom right corner
## Use spacerX and spacerY as spacers for the position of the legend, with values between 0 and 1
spacerX = 1
spacerY = 1
rect(xleft = par("usr")[2]+spacerX+1, ybottom = par("usr")[3]+spacerY+7, xright = par("usr")[2]+spacerX+1.75, ytop = par("usr")[3]+spacerY+1.5, col = "#d8d78f", lwd = 1, border = "#d8d78f")
rect(xleft = par("usr")[2]+spacerX+1.75, ybottom = par("usr")[3]+spacerY+7, xright = par("usr")[2]+spacerX+2.5, ytop = par("usr")[3]+spacerY+1.5, col = "#ceb43d", lwd = 1, border = "#ceb43d")
rect(xleft = par("usr")[2]+spacerX+2.5, ybottom = par("usr")[3]+spacerY+7, xright = par("usr")[2]+spacerX+5.17, ytop = par("usr")[3]+spacerY+1.5, col = "#c34818", lwd = 1, border = "#c34818")
rect(xleft = par("usr")[2]+spacerX+5.17, ybottom = par("usr")[3]+spacerY+7, xright = par("usr")[2]+spacerX+7.84, ytop = par("usr")[3]+spacerY+1.5, col = "#7b2300", lwd = 1, border = "#7b2300")
rect(xleft = par("usr")[2]+spacerX+7.84, ybottom = par("usr")[3]+spacerY+7, xright = par("usr")[2]+spacerX+10.5, ytop = par("usr")[3]+spacerY+1.5, col = "#421300", lwd = 1, border = "#421300")

segments(x0 = par("usr")[2]+spacerX+1, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+10.5, y1 = par("usr")[2]+spacerY+7, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+1, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+1, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+1.75, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+1.75, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+2.5, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+2.5, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+5.17, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+5.17, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+7.84, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+7.84, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)
segments(x0 = par("usr")[2]+spacerX+10.5, y0 = par("usr")[3]+spacerY+7, x1 = par("usr")[2]+spacerX+10.5, y1 = par("usr")[2]+spacerY+7.5, lwd = 1)

text(par("usr")[2]+spacerX+1, par("usr")[2]+spacerY+8.3, "0", cex = 0.8)
text(par("usr")[2]+spacerX+5.75, par("usr")[2]+spacerY+8.3, expression(italic("P")~"value"), cex = 0.8)
text(par("usr")[2]+spacerX+10.5, par("usr")[2]+spacerY+8.3, "1", cex = 0.8)

text(par("usr")[2]+spacerX+1.35, par("usr")[2]+spacerY+4.2, "0.01-0.05", cex = 0.7, srt = 90)
text(par("usr")[2]+spacerX+2.115, par("usr")[2]+spacerY+4.2, "0.05-0.10", cex = 0.7, srt = 90)
text(par("usr")[2]+spacerX+3.9, par("usr")[2]+spacerY+4.2, "0.10-0.40", cex = 0.7, srt = 90)
text(par("usr")[2]+spacerX+6.5, par("usr")[2]+spacerY+4.2, "0.40-0.70", cex = 0.7, srt = 90, col="white")
text(par("usr")[2]+spacerX+9.2, par("usr")[2]+spacerY+4.2, "0.70-1.00", cex = 0.7, srt = 90, col="white")

