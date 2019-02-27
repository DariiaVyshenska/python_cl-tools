#library(reshape)
args = commandArgs(trailingOnly=TRUE)


# cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbe-pheno-analyses/analysis/L.john-L.gasseri/inoc-expts-phenos
# Rscript meta-stats_hfhs_ref.R mic1-pheno-w-cholest.csv mic1/map.txt mic3-pheno-w-cholest.csv mic3/map.txt 


dataf1 = args[1]
mapf1 = args[2]

dataf2 = args[3]
mapf2 = args[4]


attach_map_info = function(d, m){
	data = t(read.csv(d, header = T,check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"), row.names=1))
	map = read.delim(m,header = TRUE,sep="\t",check.names=FALSE, row.names=1)
	
	#print(head(data))
	#print(head(map))

    new = cbind(map, data)
    #print(head(new))
	
	#print(levels(new$Treatment))
	return(new)
}



rel_to_ncd_mean = function(df1, df2){
		
	
	bigres = ''
	for(col in colnames(df1)[-c(1,2)]){
	    #print(col)
	    
	    res = ''
	    ncd1 = '' # expt1
	    ncd2 = '' # expt2
	    
		for (grp in c("NCD", "HFHS", "LG", "LJ")){ #levels(df1$Treatment)
		    #print(grp)
			
			sdf1 = df1[which(df1$Treatment == grp), c("Treatment",col)]
			#print(sdf1)
			
			sdf2 = df2[which(df2$Treatment == grp), c("Treatment",col)]
			#print(sdf2)
			
			if(grp == 'NCD'){
				ncd1 = mean(sdf1[, col], na.rm=T)
				ncd2 = mean(sdf2[, col], na.rm=T)
			}
			
			newvec = sdf1[, col] - ncd1
			sdf1 = cbind(sdf1, newvec)
			#print(sdf1)
			
			newvec = sdf2[, col] - ncd2
			sdf2 = cbind(sdf2, newvec)
			#print(sdf2)
			
			if (is.null(nrow(res))){
				res = rbind(sdf1, sdf2)
			} else {
				res = rbind(res, sdf1, sdf2)
			}
			
		}
		
		colnames(res) = c("Treatment", "Delete", col)
		res = res[ , c("Treatment", col)]
		#print(res)
		
		
		if(is.null(nrow(bigres))){
			bigres = res
		} else {
		    subres = res[ , col, drop=F]
			bigres = cbind(bigres, subres)
		}
	}
	
	#print(bigres)
	return(bigres)

}

rel_to_hfhs_mean = function(df1, df2){
		
	
	bigres = ''
	for(col in colnames(df1)[-c(1,2)]){
	    #print(col)
	    
	    res = ''
	    hfhs1 = '' # expt1
	    hfhs2 = '' # expt2
	    
		for (grp in c("HFHS", "NCD", "LG", "LJ")){ #levels(df1$Treatment)
		    #print(grp)
			
			sdf1 = df1[which(df1$Treatment == grp), c("Treatment",col)]
			#print(sdf1)
			
			sdf2 = df2[which(df2$Treatment == grp), c("Treatment",col)]
			#print(sdf2)
			
			if(grp == 'HFHS'){
				hfhs1 = mean(sdf1[, col], na.rm=T)
				hfhs2 = mean(sdf2[, col], na.rm=T)
			}
			
			newvec = sdf1[, col] - hfhs1
			sdf1 = cbind(sdf1, newvec)
			#print(sdf1)
			
			newvec = sdf2[, col] - hfhs2
			sdf2 = cbind(sdf2, newvec)
			#print(sdf2)
			
			if (is.null(nrow(res))){
				res = rbind(sdf1, sdf2)
			} else {
				res = rbind(res, sdf1, sdf2)
			}
			
		}
		
		colnames(res) = c("Treatment", "Delete", col)
		res = res[ , c("Treatment", col)]
		#print(res)
		
		
		if(is.null(nrow(bigres))){
			bigres = res
		} else {
		    subres = res[ , col, drop=F]
			bigres = cbind(bigres, subres)
		}
	}
	
	#print(bigres)
	return(bigres)

}



df1 = attach_map_info(dataf1, mapf1)
df2 = attach_map_info(dataf2, mapf2)

#out = rel_to_ncd_mean(df1, df2)
#write.csv(out, "phenos-mic1-mic3-rel_to_ncd_mean.csv")

out = rel_to_hfhs_mean(df1, df2)
write.csv(out, "phenos-mic1-mic3-rel_to_hfhs_mean.csv")


Expt = rep_len(1, nrow(df1))
df1 = cbind(df1, Expt)
Expt = rep_len(2, nrow(df2))
df2 = cbind(df2, Expt)
ress = rbind(df1, df2)
write.csv(ress, "phenos-mic1-mic3-merged_expts_hfhs_ref.csv")


#https://www.stat.wisc.edu/courses/st850-lindstro/handouts/blocking.pdf
#lm.out <- lm(resp ~ block + recipe,data=data)
#anova(lm.out)






# t test with default: two-sided, paired = FALSE, var.equal = FALSE,
tres = ''
for(col in colnames(out)[-1]){
	print(col)
	
	subdf = out[ , c("Treatment", col)]
	#print(subdf)
	
	hfhs = subdf[which(subdf$Treatment=="HFHS"), col]
	ncd = subdf[which(subdf$Treatment=="NCD"), col]
	lg = subdf[which(subdf$Treatment=="LG"), col]
	lj = subdf[which(subdf$Treatment=="LJ"), col]
	
	HN = t.test(hfhs, ncd)	#wilcox.test(hfhs, ncd)
	#print(HN)
	
	GH = t.test(lg, hfhs)	#wilcox.test(lg, hfhs)
	#print(GH)
	
	JH = t.test(lj, hfhs)	#wilcox.test(lj, hfhs)
	#print(JH)
	
	res = cbind(mean(hfhs), mean(ncd), mean(hfhs) - mean(ncd), HN$p.value, mean(lg), mean(lg) - mean(hfhs), GH$p.value, mean(lj), mean(lj) - mean(hfhs), JH$p.value)
	
	if(is.null(nrow(tres))){
		tres = res
	} else {
		tres = rbind(tres, res)
	}
	#print(tres)

}

rownames(tres) = colnames(out)[-1]
colnames(tres) = c("HFHS", "NCD", "HFHS - NCD", "HFHS v NCD pval", "LG", "LG - HFHS", "LG v HFHS pval", "LJ", "LJ - HFHS", "LJ v HFHS pval")
print(tres)

write.csv(tres, "phenos-mic1-mic3-results_hfhs_ref.csv")







