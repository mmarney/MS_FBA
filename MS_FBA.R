library("ggplot2")
library("readMzXmlData")
library(stringr)
library("CHNOSZ")
library("plotly")
###
# The three functions below are later used within MS_FBA.
###
add_adducts <- function(listofmasses, mode){
  # This function is used for generating adducts from the compounds list/compound_tsv file for the mode specified when running MS_FBA.
  output <- data.frame()
  if (startsWith(toupper(mode),'P') == TRUE){
    for (i in 1:length(listofmasses)){
      output[i,'M+H'] <- as.numeric(listofmasses[i]) + 1.007825
      output[i,'M+Na'] <- as.numeric(listofmasses[i]) + 22.98977
      output[i,'M+K'] <- as.numeric(listofmasses[i]) + 38.96371
      output[i,'2M+H'] <- 2 *as.numeric(listofmasses[i]) + 1.007825
      output[i,'M-H2O+H'] <- as.numeric(listofmasses[i]) - 17.00274
      output[i,'M+2H2O+H'] <-as.numeric(listofmasses[i]) + 37.02896
      output[i,'M+C2H3N+H'] <- as.numeric(listofmasses[i])+ 42.03437
      output[i,'M+C2H3N+Na'] <- as.numeric(listofmasses[i]) + 64.01632
      output[i,'M+2Na-H'] <- as.numeric(listofmasses[i]) + 44.97172
      output[i,'M+2H'] <- (as.numeric(listofmasses[i]) + 2.05165)/2
      output[i,'M+H+Na'] <- as.numeric(listofmasses[i]) + 23.9976
    }
  }
  if (startsWith(toupper(mode),'N') == TRUE){
    for (i in 1:length(listofmasses)){
      output[i,'M-H'] <- as.numeric(listofmasses[i]) - 1.007825
      output[i,'M+Na-2H'] <- as.numeric(listofmasses[i]) + 22.98977 - 2.05165
      output[i,'M+K-2H'] <- as.numeric(listofmasses[i]) + 38.96371 -2.05165
      output[i,'M+CH3COO(-)'] <-as.numeric(listofmasses[i]) + 59.013853
      output[i,'M+Cl(-)'] <- as.numeric(listofmasses[i])+ 34.968304
      output[i,'M-2H'] <- (as.numeric(listofmasses[i]) - 2.05165)/2
    }
  }
  return(output)
}
isotope_count <-  function(formula){
  # This function is used to calculate the theoretical or "expected" isotope ratio for comparison with the actual from the raw data.
  ifelse(is.na(as.numeric(str_sub(formula,1,1))),multiplier <- 1, multiplier <- as.numeric(str_sub(formula,1,1)))
  if(multiplier >1){
    formula <- str_sub(formula,2)
  }
  counts<- c(0,0,0,0,0,0,0,0)
  names(counts) <- c('C','H','O','N','S','Cl','Br','Si')
  for (i in 1:length(counts)){
    if (length(which(names(makeup(formula))== names(counts)[i]) >0)){
      counts[[i]]<-makeup(formula)[[which(names(makeup(formula))== names(counts)[i])]]
    }
  }
  #The percent natural abundance data is from the 1997 report of the IUPAC Subcommittee for Isotopic
  #Abundance Measurements by K.J.R. Rosman, P.D.P. Taylor Pure Appl. Chem. 1999, 71, 1593-1607.
  # the ratio we are multiplying our counts by is the percent isotope divide by the standard percent.
  upone <- multiplier * ((counts[['C']] *1.08)+ (counts[['H']] *0.0115) +(counts[['N']] * 0.369)+(counts[['O']] *0.038)+(counts[['S']] * 0.8006) + (counts[['Si']] * 5.07775))
  uptwo <- multiplier * ((counts[['O']] * 0.2054)+(counts[['S']] * 4.519)+ (counts[['Cl']] * 31.96094) + (counts[['Br']] * 97.2775) + (counts[['Si']] * 3.34729))
  values <- paste('M+1 :',toString(upone),'; M+2 :',toString(uptwo) , sep= '')
  return(values)
}
ranked_matches_df <- function(matches){
  # This is a function tabulates the features matched with compounds from the PyFBA model's compounds list.
  match_output = data.frame()
  for (match in 1:length(strsplit(matches,' ;; ')[[1]])){
    match_output[match,'compound'] <- str_sub(strsplit(matches,' ;; ')[[1]][match],1,str_locate(strsplit(matches,' ;; ')[[1]],' : ')[match,][['start']])
    match_output[match,'adduct']  <-  str_sub(strsplit(matches,' ;; ')[[1]][match],str_locate(strsplit(matches,' ;; ')[[1]],' : ')[match,][['end']],str_locate(strsplit(matches,' ;; ')[[1]],' ; ')[match,][['start']])
    match_output[match,'isotopic difference'] <- str_sub(strsplit(matches,' ;; ')[[1]][match],str_locate(strsplit(matches,' ;; ')[[1]],' ; ')[match,][['end']])
  }
  return(match_output)
}
###
# The function below is a convenient function that can be used to visualize the raw data MZXML files if necessary.
###
plot_mzxml_spectrum <- function(mzxml_file,index){
  plot(mzxml_file[[index]][["spectrum"]][["mass"]],mzxml_file[[index]][["spectrum"]][["intensity"]],type = 'h')
}
###
# Below is the actual MS_FBA function for data analysis.The following inputs are necessary:
# "xcms_tsvfile": From downloaded XCMS Online results, "XCMS.annotated.diffreport.tsv" file
# "files": Raw data files converted to MZXML format. Character vector of the files used to generate XCMS results
# "file_prefixes": Character vector of how sample classes were defined.
# "compounds_tsv": List of compounds from the PyFBA model. Must be formatted 'compound \t monoisotopic pass \t molecular formula \n'.
# "ion_mode": 'N' for negative mode, 'P' for positive mode. Detection mode used in the mass spec analysis.
###
MS_FBA <- function(xcms_tsvfile,files,file_prefixes,compounds_tsv,ion_mode,pvalue = 0.005,ppm =5,rtmulti=1){
  results <- read.csv(xcms_tsvfile,sep='\t')
  results <- results[which(results[,'pvalue'] <= pvalue),]
  peakgroups <- unique(results[,'pcgroup'])
  matching_df <- list() # this is going to be a list of dataframes that is the output of this function
  # this for loop separates the features detected by XCMS Online into the peak groups they are and 
  # pulls out all the crucial information from the annotated tsv 
  for (group in peakgroups){
    inter_var <- results[which(results[,'pcgroup'] == group), c('featureidx','mzmed','mzmin','mzmax','rtmed','rtmin','rtmax','maxint')]
    inter_var <- inter_var[order(inter_var[,'maxint'],decreasing = T),]
    matching_df[[paste('peakgroup',toString(group))]] <- inter_var
  }
  # after the above for loop we have the start of our dataframes. they will grow with each mzxml file that we are reading in 
  for (file in files){
    name <- strsplit(file,'.mzxml')[[1]] # the name is the name that the file is called from the instrument
    myxml <- readMzXmlFile(file)
    for(i in 1:length(matching_df)){
      inter_var <- matching_df[[i]]
      for (rownum in 1:nrow(inter_var)){
        # creating the mz window is based on the ppm value selected.
        mzmin <- inter_var[rownum,'mzmed'] - (inter_var[rownum,'mzmed'] * (ppm *0.000001)) # default. this will be an 
        mzmax <- inter_var[rownum,'mzmed'] + (inter_var[rownum,'mzmed'] * (ppm *0.000001)) # exact match since ppm = 0
        rtmin <- inter_var[rownum,'rtmin'] 
        rtmax <- inter_var[rownum,'rtmax'] 
        # creating the retention time window is dependent on the rtmulti selected.
        # the retention time difference is multipled and that time window is added on either side of 
        # the original time window
        rtdiff <- rtmax - rtmin
        rtmin <- rtmin - (rtdiff* rtmulti)
        rtmax <- rtmax + (rtdiff* rtmulti)
        # these below are set to 0 so that once something in present it will automatically be higher than 0
        # and since R starts everything at 1 even the indexs will be correct 
        highint <- 0
        highscan <- 0
        highmass_index <- 0
        for (scan in 1:length(myxml)){
          rt <- myxml[[scan]][["metaData"]][["retentionTime"]]/60 # because rt in the meta data is timed in secounds and in xcms in minitues 
          # we look for the time window first to limit the number of scans that we must evaluate
          if(rt > rtmin && rt < rtmax){
            mass_index <- which(myxml[[scan]][["spectrum"]][["mass"]] >= mzmin & myxml[[scan]][["spectrum"]][["mass"]] <= mzmax)
            if (length(mass_index) > 0 ){
              intensity <- myxml[[scan]][["spectrum"]][['intensity']][mass_index]
              if(intensity > highint){
                highint <- intensity 
                highscan <- scan 
                highmass_index <- mass_index
              }
            }
          }
        }
        isotopic_ratio <- 0
        if (highscan != 0){
          oneup_index <- which(myxml[[highscan]][['spectrum']][['mass']] - myxml[[highscan]][['spectrum']][['mass']][highmass_index] >= 1.00 
                               & myxml[[highscan]][['spectrum']][['mass']] - myxml[[highscan]][['spectrum']][['mass']][highmass_index] <=1.0066)
          if (length(oneup_index) > 0){
            isotopic_ratio <- (myxml[[highscan]][['spectrum']][['intensity']][oneup_index] / myxml[[highscan]][['spectrum']][['intensity']][highmass_index])*100
          }
        }
        # now we are adding on to the output list of dataframes and the columns will be named based on the file
        # highscan is the scan number within the mzxml file where the highest intensity for the feature mass is located 
        # highmass_index is the index in the spectrum where that feature mass is located
        # highint is the intensity of that mass, and it is the highest intensity for that feature mass
        matching_df[[i]][rownum,paste(name,'highscan')] <- highscan
        matching_df[[i]][rownum,paste(name,'highmass_index')] <- highmass_index
        matching_df[[i]][rownum,paste(name,'highint')] <- highint
        matching_df[[i]][rownum,paste(name,'isotopic_ratio')] <- isotopic_ratio
      }
    }
  }
  MLplot <- data.frame()
  LSplot <- data.frame()
  index <- 1
  for (peakgroup in 1:length(matching_df)){
    for (row in 1:nrow(matching_df[[peakgroup]])){
      # here we begin to take the average intensity of the features 
      # first we need to add the columns to the rows and set the value to zero as default
      matching_df[[peakgroup]][row,paste(file_prefixes[1],'average')] <- 0
      matching_df[[peakgroup]][row,paste(file_prefixes[2],'average')] <- 0
      matching_df[[peakgroup]][row,paste(file_prefixes[3],'average')] <- 0
      # the (m,l,s)count variables is count the number of each category are zero, 
      mcount <- length(which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[1]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) == 0 ))
      lcount <- length(which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[2]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) == 0 ))
      scount <- length(which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[3]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) == 0 ))
      # next if the counts are less then 2 (meaning 3 or more of the 5 biological replicates were present)
      # then the columns that we added will be set to the average of those intensities that were not equal to zero 
      if (mcount <= 2){
        matching_df[[peakgroup]][row,paste(file_prefixes[1],'average')] <- mean(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[1]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))[which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[1]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) != 0 )]]))
      }
      if (lcount <= 2){
        matching_df[[peakgroup]][row,paste(file_prefixes[2],'average')] <- mean(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[2]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))[which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[2]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) != 0 )]]))
      }
      if (scount <= 2){
        matching_df[[peakgroup]][row,paste(file_prefixes[3],'average')] <- mean(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[3]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))[which(as.numeric(matching_df[[peakgroup]][row,which(startsWith(colnames(matching_df[["peakgroup 2"]]),file_prefixes[3]) & endsWith(colnames(matching_df[["peakgroup 2"]]),'highint'))]) != 0 )]]))
      }
      
      MLplot[index,'M/Z'] <- matching_df[[peakgroup]][row,2]
      MLplot[index,'Intensity Difference'] <- matching_df[[peakgroup]][row,paste(file_prefixes[2],'average')] -matching_df[[peakgroup]][row,paste(file_prefixes[1],'average')]
      MLplot[index,'Peakgroup'] <- as.numeric(strsplit(names(matching_df)[peakgroup],' ')[[1]][2])
      MLplot[index,'measure points'] <- paste(file_prefixes[2],file_prefixes[1],sep = '-')
      MLplot[index,'rt'] <- matching_df[[peakgroup]][["rtmed"]][row]
      MLplot[index,'mean_isotopic_ratio'] <- mean(as.numeric(matching_df[[peakgroup]][row,which(endsWith(colnames(matching_df[["peakgroup 2"]]),'isotopic_ratio'))[-which( matching_df[[peakgroup]][row,which(endsWith(colnames(matching_df[["peakgroup 2"]]),'isotopic_ratio'))] == 0)]]))
      LSplot[index,'M/Z'] <- matching_df[[peakgroup]][row,2]
      LSplot[index,'Intensity Difference'] <- matching_df[[peakgroup]][row,paste(file_prefixes[3],'average')] -matching_df[[peakgroup]][row,paste(file_prefixes[2],'average')]
      LSplot[index,'Peakgroup'] <- as.numeric(strsplit(names(matching_df)[peakgroup],' ')[[1]][2])
      LSplot[index,'measure points'] <- paste(file_prefixes[3],file_prefixes[2],sep = '-')
      LSplot[index,'rt'] <- matching_df[[peakgroup]][["rtmed"]][row]
      LSplot[index,'mean_isotopic_ratio'] <- mean(as.numeric(matching_df[[peakgroup]][row,which(endsWith(colnames(matching_df[["peakgroup 2"]]),'isotopic_ratio'))[-which( matching_df[[peakgroup]][row,which(endsWith(colnames(matching_df[["peakgroup 2"]]),'isotopic_ratio'))] == 0)]]))
      index <- index +1 
    }
  }
  totalplot <- rbind.data.frame(MLplot,LSplot)
  totalplot[which(is.na(totalplot[,'mean_isotopic_ratio'])),'mean_isotopic_ratio'] <- 0
  # here is where the expected compounds list comes into the program 
  # Need to pull the expected compounds list from PyFBA
  compounds_list <- read.csv(compounds_tsv,sep = ',',stringsAsFactors = FALSE)
  compounds_list <- compounds_list[!duplicated(compounds_list[,1]), ]
  compound_plus_adduct <- add_adducts(compounds_list[,3],mode = ion_mode)
  rownames(compound_plus_adduct) <- compounds_list[,1]
  for (mass in unique(totalplot[,1])){
    potential <- ''
    for(index in 1:length(compound_plus_adduct[,1])){
      compound <- rownames(compound_plus_adduct)[index]
      for(adduct in colnames(compound_plus_adduct)){
        ifelse(is.na(str_locate(adduct,'2M')[1,][['start']]) == TRUE,multipler <- 1,multipler <- 2)
        if ((round(compound_plus_adduct[index,adduct],4) >= round(mass - mass * (ppm *0.000001),4)) & (round(compound_plus_adduct[index,adduct],4) <= round(mass + mass * (ppm *0.000001),4))){
          locations <- str_locate_all(str_replace(adduct,'\\(-\\)',''),c('\\+','-'))
          # these variables below are true/false questions that are counting the '+' or '-' symbols in the adduct 
          ## depending on which signs are present isotope ratio values are added or subtracted
          oneplus <- nrow(locations[[1]]) == 1
          twoplus <- nrow(locations[[1]]) == 2
          oneminus <- nrow(locations[[2]]) == 1
          twominus <- nrow(locations[[2]]) == 2
          if (twominus | twoplus == TRUE){
            ifelse(twoplus == TRUE, adduct_iso <- as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(str_replace(str_replace(adduct,'\\(-\\)',''),'-',''),'\\+',''),2)),';')[[1]][1],':')[[1]][2])
                   ,adduct_iso <- -1* (as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(str_replace(str_replace(adduct,'\\(-\\)',''),'-',''),'\\+',''),2)),';')[[1]][1],':')[[1]][2])))
          }
          if (oneplus == TRUE & oneminus == FALSE){
            adduct_iso <- as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(str_replace(adduct,'\\(-\\)',''),'\\+',''),2)),';')[[1]][1],':')[[1]][2])
          }
          if (oneminus == TRUE & oneplus == FALSE){
            adduct_iso <- -1*as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(str_replace(adduct,'\\(-\\)',''),'-',''),2)),';')[[1]][1],':')[[1]][2])
          }
          if (oneminus == TRUE & oneplus == TRUE){
            ifelse(locations[[1]][1,][['start']] < locations[[2]][1,][['start']],adduct_iso <- (as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(adduct,'\\(-\\)',''),locations[[1]][1,][['start']]+ 1, locations[[2]][1,][['start']] +1)),';')[[1]][1],':')[[1]][2])) 
                   -(as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(adduct,'\\(-\\)',''),locations[[2]][1,][['start']]+ 1)),';')[[1]][1],':')[[1]][2])),
                   adduct_iso <- -1*(as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(adduct,'\\(-\\)',''),locations[[2]][1,][['start']]+ 1, locations[[1]][1,][['start']] +1)),';')[[1]][1],':')[[1]][2])) 
                   +(as.numeric(strsplit(strsplit(isotope_count(str_sub(str_replace(adduct,'\\(-\\)',''),locations[[1]][1,][['start']]+ 1)),';')[[1]][1],':')[[1]][2])))
          }
          iso_diff <- totalplot[which(totalplot[,1] == mass),'mean_isotopic_ratio'] - (as.numeric(strsplit(strsplit(isotope_count(compounds_list[index,2]),';')[[1]][1],':')[[1]][2])*multipler) + adduct_iso
          potential <- paste(potential,compound,' : ',adduct,' ; ',iso_diff,' ;; ', sep = '')
        }
      }
    }
    # now that we have all the potential matches in one string we need to order them by the the closest
    # isotope difference to zero 
    rankings <- c()
    if (potential != ''){
      for ( i in 1:length(strsplit(potential, ' ;; ')[[1]])){
        rankings[i] <- as.numeric(str_sub(strsplit(potential, ' ;; ')[[1]][i],str_locate(strsplit(potential, ' ;; ')[[1]][i], ' ; ')[1,][['end']] +1 ,))
      }
      rankings <- order(abs(rankings))
      totalplot[which(totalplot[,1] == mass),'potential matches'] <- paste(strsplit(potential, ' ;; ')[[1]][rankings], collapse = ' ;; ')
    }
  }
  potential_matches <- totalplot[which(is.na(totalplot[,7]) == FALSE),]
  # Mplus1 and Mplus 2 are features that are potentially missannoted features that should be classified as isotopes
  # If we find a feature that falls within the isotope window and is within the same peak group as a feature that has been matched
  # we put that feature into a different data frame
  mplus1 <- c()
  for (i in 1:nrow(potential_matches)){
    mass <- potential_matches[i,1]
    mplus1 <-c(mplus1,which(totalplot[,1] -mass >= 1.000 & totalplot[,1] -mass <= 1.0066))
  }
  mplus1 <- totalplot[unique(mplus1),]
  mplus2 <- c()
  for (i in 1:nrow(potential_matches)){
    mass <- potential_matches[i,1]
    mplus2 <-c(mplus2,which(totalplot[,1] -mass >= 1.994 & totalplot[,1] -mass <= 2.008)) # still figureing out this value 
  }
  mplus2 <- totalplot[unique(mplus2),]
  no_matches <- totalplot[- c((as.integer(rownames(mplus1))),as.integer(rownames(potential_matches)),(as.integer(rownames(mplus2)))),] 
  d <- highlight_key(totalplot, ~totalplot[,'Peakgroup']) # change the second variable here to change what it will highlight 
  #             by clicking on a point.  totalplot[,1] is for m/z values but we can also do by RT or something else if we wanted 
  #             using totalplot[,'Peakgroup'] we are higlighting the points from the peakgroups.
  p <- ggplot(d,aes(x = totalplot[,'M/Z'], y = totalplot[,'Intensity Difference'], group = totalplot[,'rt'],
                    color = totalplot[,'measure points'])) + geom_point()+
    # d brings the highlight_key() information into the plot that we are making
    # aes() is how the formating of ggplots are done, within that, the first two arguments are x and y 
    # group adds the ability to show more information in the text box that appears when hovering over a point 
    # color tell the group which column to base the coloring of the points on. 
    # geom_point() tells the plot that we want dots not lines.
    labs(title = 'Intensity difference vs M/Z')+
    xlab('M/Z')+ ylab('Intensity Diference')
  # labs adds the title and xlab and y lab add the labels on the axis.
  gg <- ggplotly(p)
  # this converts from a ggplot to a plotly plot 
  highlight(gg, dynamic = FALSE, color = toRGB('Yellow'))
  # this adds the ability to click and highlight based on our highlight key info from above. 
  # cant figure out how to highlight things on a scatter 3d plot.. 
  pp <- plot_ly(totalplot, x = totalplot[,5], y = totalplot[,1])%>%
    layout(scene = list(xaxis = list(title = 'Retention Time'),
                        yaxis = list(title = 'M/Z'),
                        zaxis = list(title = 'Intensity Difference')))%>%
    add_markers(mode = 'markers', z = totalplot[,2], color = totalplot[,4],colors = c('Red','Blue'), type="scatter3d")
  # the following is the same set up as above but creates graphs for the potential matches only 
  d1 <- highlight_key(potential_matches, ~potential_matches[,'Peakgroup'])
  p1 <- ggplot(d1,aes(potential_matches[,'M/Z'], potential_matches[,'Intensity Difference'], group = potential_matches[,'rt'],
                      color = potential_matches[,'measure points'])) + geom_point()+
    xlab('M/Z')+ ylab('Intensity Diference')
  gg1 <- ggplotly(p1)
  highlight(gg1, dynamic = FALSE, color = toRGB('Yellow'))
  pp1 <- plot_ly(potential_matches, x = potential_matches[,5], y = potential_matches[,1])%>%
    layout(scene = list(xaxis = list(title = 'Retention Time'),
                        yaxis = list(title = 'M/Z'),
                        zaxis = list(title = 'Intensity Difference')))%>%
    add_markers(mode = 'markers', z = potential_matches[,2], color = potential_matches[,4],colors = c('Red','Blue'), type="scatter3d")
  # the following is the same set up as above but creates graphs for the non matches only, not including the mplus dataframes 
  d2 <- highlight_key(no_matches, ~no_matches[,'Peakgroup'])
  p2 <- ggplot(d2,aes(no_matches[,'M/Z'], no_matches[,'Intensity Difference'], group = no_matches[,'rt'],
                      color = no_matches[,'measure points'])) + geom_point()+
    xlab('M/Z')+ ylab('Intensity Diference')
  gg2 <- ggplotly(p2)
  highlight(gg2, dynamic = FALSE, color = toRGB('Yellow'))
  pp2 <- plot_ly(no_matches, x = no_matches[,5], y = no_matches[,1])%>%
    layout(scene = list(xaxis = list(title = 'Retention Time'),
                        yaxis = list(title = 'M/Z'),
                        zaxis = list(title = 'Intensity Difference')))%>%
    add_markers(mode = 'markers', z = no_matches[,2], color = no_matches[,4],colors = c('Red','Blue'), type="scatter3d")
  listofmatches <- c()
  for (i in 1:nrow(potential_matches)){
    matches <- potential_matches[i,7]
    listofmatches <- c(listofmatches,ranked_matches_df(potential_matches[i,7])[,1])
  }
  listofmatches <- unique(listofmatches)
  return(list('matching_dataframes' = matching_df, "listofmatches" = listofmatches,'totalplot' = list('potential_matches' = potential_matches, 'M_plus_one' = mplus1,'M_plus_two' = mplus2, 'no_match' = no_matches), '2D_graph' = gg, '3D_graph' = pp,'matches_2D_graph' = gg1,'matches_3D_graph' = pp1, 'no_match_2D_graph' = gg2,'no_match_3D_graph' = pp2))
}
