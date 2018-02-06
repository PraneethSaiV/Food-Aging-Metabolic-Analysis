get_values = function(dataf){
  n = rep(5,3)
  group = rep(1:3,n)

  for (i in 1:nrow(dataf)){
    da = as.numeric(as.vector(dataf[i,2:6]))
    dw = as.numeric(as.vector(dataf[i,7:11]))
    wa = as.numeric(as.vector(dataf[i,12:16]))
    
    tie = c(da,dw,wa)
    
    anova_result = anova(lm(tie ~ as.factor(group), na.action = na.omit))
    dataf$p_anova[i] = anova_result$`Pr(>F)`
    dataf$FC_dw_by_da[i] = mean(dw)/mean(da)
    dataf$FC_da_by_wa[i] = mean(da)/mean(wa)
    dataf$FC_dw_by_wa[i] = mean(dw)/mean(wa)
    aov_object = aov(tie ~ as.factor(group))
    tukey_result = TukeyHSD(aov_object)
    dataf$Tukey_dw_da[i] = tukey_result$`as.factor(group)`[10]
    dataf$Tukey_wa_da[i] = tukey_result$`as.factor(group)`[11]
    dataf$Tukey_wa_dw[i] = tukey_result$`as.factor(group)`[12]
  }
  return(dataf)
}

get_log_values = function(dataf){
  n = rep(5,3)
  group = rep(1:3,n)
  
  for (i in 1:nrow(dataf)){
    da = as.numeric(as.vector(dataf[i,2:6]))
    dw = as.numeric(as.vector(dataf[i,7:11]))
    wa = as.numeric(as.vector(dataf[i,12:16]))
    
    tie = c(da,dw,wa)
    
    anova_result = anova(lm(tie ~ as.factor(group), na.action = na.omit))
    dataf$p_anova[i] = anova_result$`Pr(>F)`
    dataf$FC_dw_by_da[i] = mean(dw) - mean(da)
    dataf$FC_da_by_wa[i] = mean(da) - mean(wa)
    dataf$FC_dw_by_wa[i] = mean(dw) - mean(wa)
    aov_object = aov(tie ~ as.factor(group))
    tukey_result = TukeyHSD(aov_object)
    dataf$Tukey_dw_da[i] = tukey_result$`as.factor(group)`[10]
    dataf$Tukey_wa_da[i] = tukey_result$`as.factor(group)`[11]
    dataf$Tukey_wa_dw[i] = tukey_result$`as.factor(group)`[12]
  }
  return(dataf)
}



cplot = function(dataframe, name){         
  png(paste0(name,'.png'), width = 1080, height = 1080)
  chart.Correlation(data.frame(dataframe))
  dev.off()
}

pcaplots = function(dataframe, name, confidence){
  dataframe = subset(dataframe,dataframe$p_anova < confidence)
  reference = gsub('@.*$','',dataframe$Compound)
  print(reference)
  selection = names(dataframe)[2:16]
  dataframe = t(dataframe[,selection])
  Principle_Components = prcomp(dataframe, center = T, scale. = T)
  Principle_Components$rotation # This
  temp_df = data.frame(Principle_Components$rotation[,1:2])
  colnames(temp_df)[1] = 'Dimension_1'
  colnames(temp_df)[2] = 'Dimension_2'
  rownames(temp_df) = reference
  temp_df = temp_df[order(abs(temp_df$Dimension_1)),][1:10,]
  png(paste(name,'_dim1.png'), width = 1200, height = 900)
  m = ggplot(temp_df, aes(x = rownames(temp_df),y = Dimension_1, fill = abs(Dimension_1)))
  m = m + geom_bar(stat = 'identity') + labs(x = 'Dimension 1', y = 'Contribution', title = 'Dimension Contributions')
  print(m)
  dev.off()
  temp_df = data.frame(Principle_Components$rotation[,1:2])
  colnames(temp_df)[1] = 'Dimension_1'
  colnames(temp_df)[2] = 'Dimension_2'
  rownames(temp_df) = reference
  temp_df = temp_df[order(abs(temp_df$Dimension_2)),][1:10,]
  png(paste(name,'_dim2.png'), width = 1200, height = 900)
  m = ggplot(temp_df, aes(x = rownames(temp_df),y = Dimension_2, fill = abs(Dimension_2))) 
  m = m + geom_bar(stat = 'identity') + labs(x = 'Dimension 2', y = 'Contribution', title = 'Dimension Contributions')
  print(m)
  dev.off()
  png(paste(name,'_scree.png'), width = 900, height = 900)
  plot(Principle_Components, type = 'l')
  dev.off()
  sample = data.frame(c('DA','DA','DA','DA','DA','DW','DW','DW','DW','DW','WA','WA','WA','WA','WA'))
  names(sample) = 'samples'
  png(paste(name,'_pca.png'), width = 900, height = 900)
  g <- ggbiplot(Principle_Components, obs.scale =0.5, var.scale = confidence, groups = as.factor(sample$samples),
                ellipse = TRUE, varname.abbrev = T, var.axes = F)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  print(g)
  dev.off()
}

#----------------------------------------------------Start Script---------------------------------------------------------#

library(PerformanceAnalytics)
library(ggbiplot)
library(ggplot2)

setwd('/media/vsppraneeth/01D3522569C0B1A0/Work/Meat Metabolics/Derico Setyabrata_Meat Metabolomics_ALL/Dry Aging Metabolomics_All')
confidence = 0.05

raw_data = read.csv('Dry Aging Metabolomics_ ALL (no 12)_Raw.tsv', sep = '\t')
norm_data = read.csv('Dry Aging Metabolomics_ALL (no 12)_Normalized.tsv', sep = '\t')
log_data = raw_data
log_data[,2:16] = log(log_data[,2:16],2)


# Obtain Values
raw_data = get_values(raw_data)
norm_data = get_values(norm_data)
log_data = get_log_values(log_data)

# Correlation Plots
# Raw
cplot(raw_data[,2:6],'Raw_DA')
cplot(raw_data[,7:11],'Raw_DW')
cplot(raw_data[,12:16],'Raw_WA')

# Norm
cplot(norm_data[,2:6],'Norm_DA')
cplot(norm_data[,7:11],'Norm_DW')
cplot(norm_data[,12:16],'Norm_WA')

# Log
cplot(log_data[,2:6],'Log_DA')
cplot(log_data[,7:11],'Log_DW')
cplot(log_data[,12:16],'Log_WA')

# PCA Plots (For only significant)
pcaplots(raw_data,'Raw', 0.05)
pcaplots(norm_data,'Norm', 0.05)
pcaplots(log_data,'Log', 0.05)


# Save Files
write.table(raw_data, 'Raw_Data.tsv', sep = '\t', row.names = FALSE)
write.table(norm_data, 'Norm_Data.tsv', sep = '\t', row.names = FALSE)
write.table(log_data, 'Log_Data.tsv', sep = '\t', row.names = FALSE)











