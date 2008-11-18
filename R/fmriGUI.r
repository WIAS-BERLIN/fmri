fmri.gui <- function() {
        require(tcltk) || stop("install package tcltk!")
	
	# functions

	# select file, name saved in tclVar dataFile
	selectFirstDataFile <- function(){
        	tclvalue(dataFileFirst) <- tcl("tk_getOpenFile")  
	}	
	selectDesignFileSave <- function(){
        	tclvalue(designFile) <- tcl("tk_getSaveFile")  
	}
	selectDesignFileLoad <- function(){
        	tclvalue(designFile) <- tcl("tk_getOpenFile")
	}

	# view p values in 3D
	results3D <- function(){
		tkmessageBox(message="Good Choice! But still under construction ... ")
	}
	
	# view p values as slices
	resultsSlices <- function(){
		ttt <- extract.data(data)
		plot(pvalue,anatomic=ttt[,,,1])
	}

	# save design in filepath/filename.Design
	# format: (number of scans, sum of number of onsets for all conditions, number of onsets (n1) for condition1 (of M), first onset condition1, second onset, ..., 	
        #          nM-th onset conditionM, sum of number of durations for all conditions, number of durations (d1) for condtion1, first duration condition1, ..., dM-th duration 
	#	   conditionM, interscan intervals)  
	saveDesign <- function(){
		print("Saving ...")		
		onsVec <- c(0)
		durVec <- c(0)
		for (i in 1:length(onsets)){ # number of all onsets
			onsVec[i] = onsets[i]
		}
		for (i in 1:length(durations)){ # number of all explicitely given durations
			durVec[i] = durations[i]
		}			
		textVar <- c(0)
		gap = "" # gap between entries, otherwise chaos ;)
		textVar[1] = scanspS # number of scans
		textVar[2] = gap
		textVar[3] = length(onsets)
		textVar[4] = gap	
		for (i in 1:length(onsets)){
			textVar[2*i+3] = onsets[[i]] # onsets
			textVar[2*i+4] = gap	
		} 
		currPos = 2*length(onsets)+4
		textVar[currPos+1] = length(durations)
		textVar[currPos+2] = gap
		currPos = currPos+2
		for (i in 1:length(durations)){
			textVar[2*i+currPos-1] = durations[[i]] # durations
			textVar[2*i+currPos] = gap
		} 
		currPos = currPos + 2*length(durations)
		textVar[currPos+1] = interscanInt # interscan intervals
		textVar[currPos+2] = gap	
		textVar[currPos+3] = as.character(tclvalue(rbDesign))
		designFileOld = as.character(tclvalue(designFile))
		selectDesignFileSave()
		if (tclvalue(designFile) != ""){
			designFileText <- as.character(tclvalue(designFile))
			write(textVar,designFileText) # design written	
			print("Saving completed")
		}
	}

	loadDesignHelp <- function(){
		if (nrStep>=1){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startDataLayout(loadPlan=1)
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to reload the design? \n All progresses will be lost.",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startDataLayout(loadPlan=1)	
	}
	
	createDesignHelp <- function(){
		if (nrStep>=1){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startDataLayout(loadPlan=0)
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to recreate the design? \n All progresses will be lost.",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startDataLayout(loadPlan=0)	
	}

	# load Design, saved in filepath/dataFile.Design
	# format (see above) deciphered and save partly as reals, partly as real-matrices (durations/onsets of condition i in column i saved, and some placeholders (-1))	
	loadDesign <- function(){
		print("Loading ...")

		selectDesignFileLoad()
			
		if (tclvalue(designFile)==""){
			nrStep <<- nrStep - 1
		}
		else {
			designFileText <- as.character(tclvalue(designFile))
			text = readChar(designFileText,1000) # read file
			designData = unlist(strsplit(text,"\n\n")) # delete "\n\n"
			designData[length(designData)] = strsplit(designData[length(designData)],"\n") # delete last "\n"
			scanspS <<- as.numeric(designData[1]) # number of scans
			interscanInt <<- as.numeric(designData[length(designData)-1]) # interscan intervals
			tclvalue(rbDesign) <<- as.character(designData[length(designData)])
			lengthOns = as.numeric(designData[2])
			lengthDur = as.numeric(designData[2+lengthOns+1])
			i=0
			loops = 0
			onsetsMat = matrix(-1,lengthOns,lengthOns)	
			textOns = list()
			textDur = list()	
			while (i  < lengthOns){
				textCurr = ""
				loops = loops + 1
				lengthCurr = as.numeric(designData[3+i])
				onsets[i+1] <<- lengthCurr
				for (j in 1:lengthCurr){
					onsets[i+j+1] <<- as.numeric(designData[3+i+j]) # onsets as vector
					textCurr = paste(textCurr,"",onsets[i+j+1])
				}
				i = i + lengthCurr + 1
				textOns[loops] = textCurr
				nrCond <<- loops
			}
			i=0
			loops = 0
			durationsMat = matrix(-1,lengthDur,lengthDur)
			while (i  < lengthDur){
				textCurr = ""
				loops = loops + 1
				lengthCurr = as.numeric(designData[3+i+lengthOns+1])
				durations[i+1] <<- lengthCurr
				for (j in 1:lengthCurr){
					durations[i+j+1] <<- as.numeric(designData[3+i+j+lengthOns+1]) # durations as matrix
					textCurr = paste(textCurr,"",durations[i+j+1])
				}
				i = i + lengthCurr + 1
				textDur[loops] = textCurr
			}
			tclvalue(rbDesign) <<- as.character(designData[length(designData)])
			print("Loading completed")
			#Print read in design for control		
			text1 = paste("Interscan intervals",":",interscanInt)
			text2 = paste("Scans per session",":",scanspS)
			text3 = paste("Number of conditions",":",nrCond)
			text4 = paste("Design in",":",as.character(tclvalue(rbDesign)))	
			print(text1)
			print(text2)
			print(text3)
			print(text3)
			for (i in 1:nrCond){
				textC = paste("Condition","",i)
				textO = paste("Onsets",":",textOns[i])
				textD = paste("Durations",":",textDur[i])
				textTog = paste(textO," - ",textD)
				print(textC)
				print(textTog)
			}	
			startData()	
		}
	}

	# design created
	# interscan intervals, number of scans, number of conditions, onsets of the conditions and durations defined
	# function started from the start window (base.aws)
	createDesign <- function(){	
		# new window
		if(.Platform$OS.type == "windows") flush.console()
    		base.design <- tktoplevel(bg=wiasblue)
    		tkwm.title(base.design, "FMRI Analysis - Design Definition")
	
		# functions
			
		# quits design window
		design.quit <- function(...){
			onYes <- function(){
				tkdestroy(base.design)
				tkdestroy(ttWarning)
				nrStep <<- nrStep-1
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to quit the design definition? ",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)	
    		}
	
		startDataHelp <- function(){
			tkdestroy(base.design)
			startData()
		}
		
		# condition attributes defined here
		# called from design.create
		design.conditions <- function(){
			tkdestroy(frame4)
			n.comp = as.numeric(tclvalue(nrCondTc))
			globali = 0			
		
			# functions
				
			# creates new line with entries for condition's name, onsets and durations			
			nextCondition <- function(){
				#assign, newenv
				globali <<- globali+1
				# conditon data read in
				readCond <- function(){					
					durationsTc[globali] <<- tclvalue(currDur)
					onsetsTc[globali] <<- tclvalue(currOnset)	
					namesTc[globali] <<- tclvalue(currName)			
					nextCondition()
				}
				# called for the last condition
				# data read in and button for start smoothing added to condition window (calls function smoothing())
				condition.next <- function(){
					durationsTc[globali] <<- tclvalue(currDur)
					onsetsTc[globali] <<- tclvalue(currOnset)
					namesTc[globali] <<- tclvalue(currName)	
					tkgrid(OkButton,QuitButton,SaveButton,padx=10,pady=10)	
					# umformungen !!!!	
					#Interscan intervals, number of scans and number of conditions saved as numerics		
					interscanInt <<- as.numeric(tclvalue(interscanIntTc))
					scanspS <<- as.numeric(tclvalue(scanspSTc))
					nrCond <<- as.numeric(tclvalue(nrCondTc))	
					nrCondCurr <- nrCond
					# help vars				
					onsetsHelp <- list()
					durationsHelp <- list()
					onsSpl = list()
					durSpl = list()
					globali = 1
					globalj = 1
					onsVecCurr2 <- c(0)
					durVecCurr2 <- c(0)					
					onsVecCurr <- c(0)
					durVecCurr <- c(0)
					onsSpl2 = list()
					durSpl2  = list()
					# writing onsets and durations in respectively a vector (as numerics), using the following format
					# (condition1 number of onsets (n), first onset, second onset, ..., n-th onset, condition2 number of onsets, ...)
					# the same with durations
			
					for (i in 1: n.comp){
						onsetsHelp[[i]] <- as.vector(onsetsTc[i]) #getting the vectors (char)
						durationsHelp[[i]] <- as.vector(durationsTc[i])
						# for each vector of onsets and durations separator determined (" ","," or ";")
						sepO = " "
						allCharsO = unlist(strsplit(onsetsHelp[[i]],""))
						print(onsetsHelp)							
						print(allCharsO)
						for (j in 1:length(allCharsO)){
							if (allCharsO[j]==","){
								sepO = ","
							}
							if (allCharsO[j]==";"){
								sepO = ";"					
							}
						}
						sepD = " "
						allCharsD = unlist(strsplit(durationsHelp[[i]],""))
						for (j in 1:length(allCharsD)){
							if (allCharsD[j]==","){
								sepD = ","
							}
							if (allCharsD[j]==";"){
								sepD = ";"					
							}
						}
						onsSpl[i] <- strsplit(onsetsHelp[[i]],sepO) # split with above found separator	
						durSpl[i] <- strsplit(durationsHelp[[i]],sepD) # split with above found separator
						#currInd <<- 1						
						#for (j in 1:length(onsSpl2[i])){
						#	if (onsSpl2[i,j]!=""){ 	
						#		print(currInd)	
						#		onsSpl[i,currInd] <<- onsSpl2[i,j] 
						#		currInd <<- currInd + 1
						#	}	
						#}		
						onsVecCurr2 <- unlist(onsSpl[i]) # unlist entries
						durVecCurr2 <- unlist(durSpl[i])
						
						currIndO <<- 1	
						for (j in 1:length(onsVecCurr2)){
							if (onsVecCurr2[j]!=""){ 	
								print(currIndO)	
								onsVecCurr[currIndO] <- onsVecCurr2[j] 
								currIndO <<- currIndO + 1
							}	
						}
						
						currIndD <<- 1	
						for (j in 1:length(durVecCurr2)){
							if (durVecCurr2[j]!=""){ 	
								print(currIndD)	
								durVecCurr[currIndD] <- durVecCurr2[j] 
								currIndD <<- currIndD + 1
							}	
						}
							
						print(onsSpl[i])
						print(onsVecCurr)
						onsets[globali] <<- length(onsVecCurr)
						durations[globalj] <<- length(durVecCurr)
						for (j in 1:length(onsVecCurr)){
							onsets[j+globali] <<- as.numeric(onsVecCurr[j]) #single entries transformed to numerics and saved
						} 
						for (j in 1:length(durVecCurr)){
							durations[j+globalj] <<- as.numeric(durVecCurr[j])
						} 	
						globali <- globali + 1 + length(onsVecCurr)	
						globalj <- globalj + 1 + length(durVecCurr)
					}			
				}	
				currOnset <- tclVar()
				currDur <- tclVar()
				currName <- tclVar("Condition Name")
				labelText = paste("Condition","",globali)
				condLabel  <- tklabel(tkCurr,text = labelText,bg=wiasblue,width=20,font="Arial 12 bold")	
				nameLabel <- tklabel(tkCurr,text="Name",bg=wiasblue,font="Arial 12 bold")	
				onsetsLabel <- tklabel(tkCurr,text="Onset times",bg=wiasblue,font="Arial 12 bold")	
				durationsLabel <- tklabel(tkCurr,text="Duration",bg=wiasblue,font="Arial 12 bold")
				nameEntry = tkentry(tkCurr,textvariable=currName,width=25,bg="#ffffff")
				onsetsEntry = tkentry(tkCurr,textvariable=currOnset,width=25,bg="#ffffff")
				durationsEntry = tkentry(tkCurr,textvariable=currDur,width=25,bg="#ffffff")
				cOkB <- tkbutton(tkCurr,text="Next Condition",width=15,command=readCond,bg=wiaslightblue)
				cReadyB <- tkbutton(tkCurr,text="Ready",width=15,command=condition.next,bg=wiaslightblue)		
				if (globali < n.comp){
					tkgrid(condLabel,nameLabel,nameEntry,padx=10,pady=10)
					tkgrid(onsetsLabel,onsetsEntry,padx=10,pady=10)	
					tkgrid(durationsLabel,durationsEntry,cOkB,padx=10,pady=10)
					tkgrid.configure(onsetsLabel,column=1)
					tkgrid.configure(onsetsEntry,column=2)
					tkgrid.configure(durationsLabel,column=1)
					tkgrid.configure(durationsEntry,column=2)
					tkgrid.configure(cOkB,column=3)
				}
				if (globali == n.comp){
					tkgrid(condLabel,nameLabel,nameEntry,padx=10,pady=10)
					tkgrid(onsetsLabel,onsetsEntry,padx=10,pady=10)	
					tkgrid(durationsLabel,durationsEntry,cReadyB,padx=10,pady=10)
					tkgrid.configure(onsetsLabel,column=1)
					tkgrid.configure(onsetsEntry,column=2)
					tkgrid.configure(durationsLabel,column=1)
					tkgrid.configure(durationsEntry,column=2)
					tkgrid.configure(cReadyB,column=3)	
				}				
			}
	
			# condition window created
			tkCurr = tkframe(base.design,relief="groove",borderwidth=2,bg=wiasblue)		
				
			nextCondition()	# next condition started -> first line for condition attributes created 
			
			tkgrid(tkCurr,sticky="ew")
			tkgrid(frame5)	
		}		

		# frames in the design window (base.design)
	
		# frame for defining the time between scans and the number of scans
		frame1D <- tkframe(base.design,relief="groove",borderwidth=2,bg=wiasblue)
		objIntervalL <- tklabel(frame1D,text="Interscan Intervals",bg=wiasblue,font="Arial 12 bold")	
		objIntervalE <-	tkentry(frame1D,textvariable=interscanIntTc,width=6,bg="#ffffff")
		tkgrid(objIntervalL,objIntervalE,padx=10,pady=10)
		objscanspSTcL <- tklabel(frame1D,text="Scans per Session",bg=wiasblue,font="Arial 12 bold")	
		objscanspSTcE <- tkentry(frame1D,textvariable=scanspSTc,width=6,bg="#ffffff")
		tkgrid(objscanspSTcL,objscanspSTcE,padx=10,pady=10)
		tkgrid(frame1D,sticky="ew")
		
		# frame for the choice of designing in scans or seconds
		frame2D <- tkframe(base.design,relief="groove",borderwidth=2,bg=wiasblue)
		rbD1 <- tkradiobutton(frame2D,bg=wiasblue)
		rbD2 <- tkradiobutton(frame2D,bg=wiasblue)
		tkconfigure(rbD1,variable=rbDesign,value="scans",bg=wiasblue)
		tkconfigure(rbD2,variable=rbDesign,value="seconds",bg=wiasblue)
		tkgrid(tklabel(frame2D,text="Design in scans or seconds?",bg = wiasblue,font="Arial 12 bold"))
		tkgrid(tklabel(frame2D,text="scans ",bg = wiaslightblue,font="Arial 12"),rbD1,tklabel(frame2D,text="seconds ",bg = wiaslightblue,font="Arial 12"), rbD2, padx=10,pady=10)
		tkgrid(frame2D,sticky="ew")	
	
		# frame for appointing the number of conditions
		frame3 <- tkframe(base.design,relief="groove",borderwidth=2,bg=wiasblue)
		objNrCondL <- tklabel(frame3,text="Number of conditions",bg=wiasblue,font="Arial 12 bold")	
		objNrCondE <- tkentry(frame3,textvariable=nrCondTc,width=6,bg="#ffffff")
		tkgrid(objNrCondL,objNrCondE,padx=10,pady=10)
		tkgrid(frame3,sticky="ew")
	  	
		# frame for "Ok" and Quit" button
		# pushing the "Ok" button will start design.condition (see above)
		# frame quitted after defining all conditions
    		frame4 <- tkframe(base.design, borderwidth = 2, bg = wiasblue)
    		ok.but<- tkbutton(frame4, text = "Ok", bg = wiaslightblue, command = design.conditions)
    		q.but <- tkbutton(frame4, text = "Quit", command = design.quit, bg = wiaslightblue)
    		tkgrid(ok.but, q.but, padx = 30, pady = 20)
    		tkgrid(frame4)
	
		# frame for the option to save the design
		frame5 <- tkframe(base.design,relief="groove",borderwidth=2,bg=wiasblue)
		OkButton <- tkbutton(frame5,text="Ok",font="Arial 12",command=startDataHelp,bg = wiaslightblue)
		QuitButton <- tkbutton(frame5, text = "Quit",font="Arial 12", command = design.quit, bg = wiaslightblue)
		SaveButton <-tkbutton(frame5, text = "Save",font="Arial 12", bg = wiaslightblue, command = saveDesign)
			
	}

	startDataLayout <- function(loadPlan=-1){
		onsets <<- list()
		durations <<- list()
		layoutFunction(1)
		if (loadPlan==0) createDesign()
		if (loadPlan==1) loadDesign()
	}

	startData <- function(){	
		tkgrid(objFileL,pady = 10, padx = 10,sticky="ew")
		tkgrid(frameData1)
		tkgrid(objFileB1,objFileE1,pady = 10, padx = 10)
		tkgrid(frameData2)
		tkgrid(objFileLoad,helpButton2,pady=10,padx=15)
		tkgrid(frameData4)
		tkgrid(helpLabel1)	
		tkgrid(frameData3,sticky="ew")		
	}

	startContrastHelp <- function(){
		if (nrStep>=2){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startContrast()
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to readjust the mask? \n All subsequent steps have to be redone", font="Arial 13", bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startContrast()	

	}

	startContrast <- function(){
		layoutFunction(2)	
		# Read in data(AFNI, ANALYZE or NIFTI file)
		# Report about successful read in

		quantile = as.numeric(tclvalue(quantileTc))
		data$mask=extract.data(data)[,,,1]>quantile
	
		tkgrid(objcontrastL,padx=10,pady=10)
		tkgrid(frameContrast1)	
		tkconfigure(objcontrastE,state="normal")
		tkgrid(objcontrastE,objcontrastB,helpButton4,padx=10,pady=10)
		tkgrid(frameContrast2)	
		
		if (nrCond==1){
			tkconfigure(objcontrastE,state="readonly")
			startEstimationHelp()			
		}
		
	}

	startEstimationHelp <- function(){
		if (nrStep>=3){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startEstimation()
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to change the contrast? \n All subsequent steps have to be redone.",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startEstimation()	
	}

	startEstimation <- function(){
		layoutFunction(3)
		
		tkgrid(helpLabel2)
		tkgrid(frameContrast3,stick="ew")
		# main evaluations:
		# 1.) data <- read.AFNI(tclvalue(dataFile))	
		# 2.) hrfCurr = fmri.stimulus(nrSc,onsCurrCond,durCurrCond,intSc)
		# 3.) x <- fmri.design(hrf)
		# 4.) spm <- fmri.lm(data,x)
		# 5.-7.) cf. function smooth	
		#
		mycontrast <<- c()

		if (nrCond==1){
			tclvalue(mycontrastTc) = 1
			mycontrast[1]<<-1
		}
		else {			
			mycontrastHelp <- list()
			contrSplVec <- c(0)	
	
			mycontrastHelp <- as.vector(tclvalue(mycontrastTc))
			allCharsC = unlist(strsplit(mycontrastHelp,""))
			sepC = " "		
			for (j in 1:length(allCharsC)){
				if (allCharsC[j]==","){
					sepC = ","
				}
				if (allCharsC[j]==";"){
					sepC = ";"					
				}
			}	
			contrSplVec <- unlist(strsplit(mycontrastHelp,sepC)) 
			
			for (j in 1:length(contrSplVec)){
				mycontrast[j] <<- as.numeric(contrSplVec[j]) # single entries transformed to numerics and saved
			} 
		}

		# create expected bold signal for each condition
		hrf = list()
		globalpos1 = 1
		globalpos2 = 1
		for (i in 1:nrCond){
			onsCurrCond = c(0)	
			durCurrCond = c(0)
			for (j in 1:(onsets[[globalpos1]])){
				onsCurrCond[j] = as.numeric(onsets[[globalpos1+j]])
			}
			for (j in 1:durations[[globalpos2]]){
				durCurrCond[j] = as.numeric(durations[[globalpos2+j]])
			}

			print(durCurrCond)
			print(onsCurrCond)

			if (i==1){
				if (as.character(tclvalue(rbDesign))=="scans"){
					hrf = fmri.stimulus(scanspS,onsCurrCond,durCurrCond,interscanInt)			
				}
				else { # secs
					hrf = fmri.stimulus(scanspS,times=onsCurrCond,durCurrCond,interscanInt)
				}
			}
			else {
				if (as.character(tclvalue(rbDesign))=="scans"){
					hrfCurr = fmri.stimulus(scanspS,onsCurrCond,durCurrCond,interscanInt)
				}
				else {	# secs
					hrfCurr = fmri.stimulus(scanspS,times=onsCurrCond,durCurrCond,interscanInt)
				}			
				hrf = cbind(hrf,hrfCurr)
			}				
			globalpos1 = globalpos1 + onsets[[globalpos1]] + 1	
			globalpos2 = globalpos2 + durations[[globalpos2]] + 1	
		}

		# create design matrix
		# report about successful creation 
		x <- fmri.design(hrf)	
		print("Stimulus calculated")		
			
		# generate parametric map from linear model
		# report about successful generation 	
		spm <<- fmri.lm(data,x,contrast=mycontrast)
		print("Parametric map calculated")		
				
		
		delta = data$delta
		tclvalue(hMax) = round(100*as.numeric(10/min(delta)))/100

		tkgrid(smoothChoiceL,padx=10,pady=10)
		tkgrid(frameSmoothChoice1)
		tkgrid(objhMaxL,objhMaxE,objButtonSmooth,padx=13,pady=10)
		tkgrid(frameSmoothChoice2)
		tkgrid(smoothChoiceB2,helpButton5,padx=15,pady=10)
		tkgrid(frameSmoothChoice4)
		tkgrid(helpLabel4)
		tkgrid(frameSmoothChoice3,sticky="ew")		
	}
		
	contWithoutHelp <- function(){
		if (nrStep>=4){		
			onYes <- function(){
				tkdestroy(ttWarning)
				contWithout()
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to continue? \n The p values will get recalculated.",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else contWithout()	
	}

	contWithout <- function(){
		layoutFunction(4)
		# calculation of the p-values and reporting about its success
		pvalue <<-fmri.pvalue(spm)
		print("p values calculated")
	
		# choice to view the p-values in 3d or as slices
		tkgrid(objResultsL,objResultsB1,objResultsB2,helpButton6,padx=10,pady=10)
		tkgrid(frameResults,sticky="ew");	
	}

	startSmoothingHelp <- function(){
		if (nrStep>=4){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startSmoothing()
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to (re)smooth the parametric map? \n The p values will get recalculated, too.", font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startSmoothing()	
	}

	startSmoothing <- function(){
		layoutFunction(4)
		print(as.numeric(tclvalue(hMax)))
			
		# smoothing of the parametric map and reporting about its success
		spmsmooth <- fmri.smooth(spm,hmax=as.numeric(tclvalue(hMax)))
					
		print("Parametric map adaptive smoothed")
	
	
		# calculation of the p-values and reporting about its success
		pvalue <<-fmri.pvalue(spmsmooth)
		print("P values calculated")	
		
		# choice to view the p-values in 3d or as slices
		tkgrid(objResultsL,objResultsB1,objResultsB2,padx=10,pady=10)
		tkgrid(frameResults,sticky="ew");				
	}

	layoutFunction <- function(currStep){
		if (nrStep>=1 && currStep<=1){	
			tkgrid.forget(frameData1)
			tkgrid.forget(frameData2)
			tkgrid.forget(frameData4)
			tkgrid.forget(frameData3)
			tkgrid.forget(helpLabel1)
		}	
		
		if (nrStep>=1.5 && currStep<=1.5){
			if (conform == 1){				
				tkgrid.forget(frameMask1)
				tkgrid.forget(frameMask2)
				tkgrid.forget(frameMask4)
				tkgrid.forget(helpLabel5)
				tkgrid.forget(frameMask3)			
			}		
		}
		if (nrStep>=2 && currStep<=2){				
			tkgrid.forget(frameContrast1)	
			tkgrid.forget(frameContrast2)
			tkgrid.forget(helpLabel2)
			tkgrid.forget(frameContrast3)	
		}
		if (nrStep>=3 && currStep<=3){
			tkgrid.forget(frameSmoothChoice1)
			tkgrid.forget(frameSmoothChoice2)
			tkgrid.forget(frameSmoothChoice4)
			tkgrid.forget(helpLabel4)
			tkgrid.forget(frameSmoothChoice3)
		}
		if (nrStep>=4 && currStep<=4){
			tkgrid.forget(frameResults)			
		}

		nrStep <<- currStep		
	}	

	viewMask <- function(width=14,height=7){
		ttt <- extract.data(data)
		dim <- data$dim		
		X11(width=12,height=7)
		nrrow <<- 0
		nrcol <<- 0
		if (round(sqrt(dim[3]))==sqrt(dim[3])){
			nrrow <<- sqrt(dim[3])
			nrcol <<- sqrt(dim[3])	
		}
		else {
			if ((ceiling(sqrt(dim[3]))-1)*ceiling(sqrt(dim[3])) >= dim[3]) {
				nrrow <<- ceiling(sqrt(dim[3]))-1
				nrcol <<- ceiling(sqrt(dim[3]))
			}
			else {
				nrrow <<- ceiling(sqrt(dim[3]))
				nrcol <<- ceiling(sqrt(dim[3]))
				
			}
		}
		mat = matrix(0,nrrow,nrcol+1)
		for (i in 1:nrrow){ for (j in 1:(nrcol+1)){ if ((i-1)*(nrcol+1)+j-(i-1) <= dim[3]) { mat[i,j]=(i-1)*(nrcol+1)+j-(i-1) } } }
		for (i in 1:nrrow){ mat[i,nrcol+1] = dim[3]+1 }
		print(mat)
		widthsvec = c(1:nrcol+1)
		for (i in 1:nrcol){ widthsvec[i]=0.5/nrcol }
		widthsvec[nrcol+1] = 0.5
		layout(mat,widthsvec)
		par(mar=c(0.5,0.5,0.5,0.5))
		for (i in 1:dim[3]) image(ttt[,,i,1]>as.numeric(tclvalue(quantileTc)),yaxt="n",xaxt="n")
		par(mar=c(5,5,3,1))
		d0 <- density(ttt[,,,1])	
		d1 <- density(ttt[round((1/8)*dim[1]):round((7/8)*dim[1]),round((1/8)*dim[2]):round((7/8)*dim[2]),round((1/8)*dim[3]):round((7/8)*dim[3]),1])
		d2 <- density(ttt[round((2/8)*dim[1]):round((6/8)*dim[1]),round((2/8)*dim[2]):round((6/8)*dim[2]),round((2/8)*dim[3]):round((6/8)*dim[3]),1])	
		d3 <- density(ttt[round((3/8)*dim[1]):round((5/8)*dim[1]),round((3/8)*dim[2]):round((5/8)*dim[2]),round((3/8)*dim[3]):round((5/8)*dim[3]),1])	
		plot(d0)
		lines(d1,col=2)
		lines(d2,col=3)
		lines(d3,col=4)		
		lines(c(as.numeric(tclvalue(quantileTc)),as.numeric(tclvalue(quantileTc))),range(d0$y),col=6)
		maxi <<- c(1:4)
		maxVal = c(1:4)
		maxVal[4] = max(d0$y)
		maxVal[1] = max(d1$y)
		maxVal[2] = max(d2$y)
		maxVal[3] = max(d3$y)
		for (i in 1:length(d0$y)){ if (d0$y[i]==maxVal[4]){ maxi[4] <<- i } }
		for (i in 1:length(d1$y)){ if (d1$y[i]==maxVal[1]){ maxi[1] <<- i } }
		for (i in 1:length(d2$y)){ if (d2$y[i]==maxVal[2]){ maxi[2] <<- i } }
		for (i in 1:length(d3$y)){ if (d3$y[i]==maxVal[3]){ maxi[3] <<- i } }
		text4 = paste("Density of data")			
		text1 = paste("Density of data[",round((1/8)*dim[1]),":",round((7/8)*dim[1]),",",round((1/8)*dim[2]),":",round((7/8)*dim[2]),",",sep="")
		text1 = paste(text1,round((1/8)*dim[3]),":",round((7/8)*dim[3]),"]",sep="")		
		text2 = paste("Density of data[",round((2/8)*dim[1]),":",round((6/8)*dim[1]),",",round((2/8)*dim[2]),":",round((6/8)*dim[2]),",",sep="")
		text2 = paste(text2,round((2/8)*dim[3]),":",round((6/8)*dim[3]),"]",sep="")
		text3 = paste("Density of data[",round((3/8)*dim[1]),":",round((5/8)*dim[1]),",",round((3/8)*dim[2]),":",round((5/8)*dim[2]),",",sep="")
		text3 = paste(text3,round((3/8)*dim[3]),":",round((5/8)*dim[3]),"]",sep="")
		text(d1$x[maxi[1]],d1$y[maxi[1]],text1,adj=c(0,0),col=2,cex=1.5)
		text(d2$x[maxi[2]],d2$y[maxi[2]],text2,adj=c(0,0),col=3,cex=1.5)		
		text(d3$x[maxi[3]],d3$y[maxi[3]],text3,adj=c(0,0),col=4,cex=1.5)
		text(d0$x[maxi[4]],0.8*d0$y[maxi[4]],text4,adj=c(0,0),col=1,cex=1.5)
		text(as.numeric(tclvalue(quantileTc)),0.5*d0$y[maxi[4]],"Threshold",adj=c(0,0),col=6,cex=1.5)
	}

	startAdjustHelp <- function(){
		if (nrStep>=1.5){		
			onYes <- function(){
				tkdestroy(ttWarning)
				startAdjust()
			}	
			onNo <- function(){
				tkdestroy(ttWarning)
			}
			ttWarning = tktoplevel(bg=wiasblue)
			tkwm.title(ttWarning, "Affirmation")
			warningFrame1 = tkframe(ttWarning,bg=wiasblue)
			warningLabel = tklabel(warningFrame1,text="Are you sure you want to reload the data? \n All progresses (out of design definitions) will be lost.",font="Arial 13",bg=wiasblue)
			warningFrame2 = tkframe(ttWarning,bg=wiasblue)	
			warningB1 = tkbutton(warningFrame2,text="Yes",command=onYes)
			warningB2 = tkbutton(warningFrame2,text="No",command=onNo)	
			tkgrid(warningLabel)
			tkgrid(warningB1,warningB2,padx=10,pady=10)
			tkgrid(warningFrame1)	
			tkgrid(warningFrame2)
		}  else startAdjust()
	}

	startAdjust <- function(){
		layoutFunction(1.5)

		dataFile = as.character(tclvalue(dataFileFirst))
		help <- tolower(unlist(strsplit(dataFile,"")))
		help2 <- unlist(strsplit(dataFile,""))
		nrChars <- length(help) # nchar(dataFile)		
		if (help[nrChars-2]=="i"&& help[nrChars-1]=="m" && help[nrChars]=="g"){
			dataType = "ANALYZE"
		}				
		else if (help[nrChars-2]=="n"&& help[nrChars-1]=="i" && help[nrChars]=="i"){
			dataType = "NIFTI"
		}
		else if ((help[nrChars-3]=="h"&& help[nrChars-2]=="e" && help[nrChars-1]=="a" && help[nrChars]=="d")||(help[nrChars-3]=="b"&& help[nrChars-2]=="r" && help[nrChars-1]=="i" && help[nrChars]=="k")){
			dataType = "AFNI"
		}
		else if (help[nrChars-2]=="h"&& help[nrChars-1]=="d" && help[nrChars]=="r"){
			if (file.info(dataFile)$size == 348) dataType = "ANALYZE"
			else dataType = "NIFTI"
		}	
		else {
			dataType = "unknown"	
		}
		
		dataTypeGlobal = dataType			

		if (dataType=="AFNI"){
			data <<- read.AFNI(as.character(tclvalue(dataFileFirst)))
		}

		if (dataType=="ANALYZE"){
			index <<- -1
			indexhelp <<- -1	
			for (helpi in 1: nrChars-2){
				i = nrChars-1-helpi	
				if (index==-1){
					isFigure <<- FALSE	
					for (j in 0:9){
						if (help[i]==as.character(j)) {
							isFigure <<- TRUE
						}				
					}
					if (isFigure && indexhelp==-1) indexhelp <<- i
					if (!isFigure && indexhelp!=-1) indexhelp <<- -1
					if (isFigure && indexhelp == i+2) index <<- indexhelp-2
				}
			}
			if (index!=-1){
				myprefix = ""
				mypostfix= ""
				for (i in i:index-1){
					myprefix = paste(myprefix,help2[i],sep="")
				}
				lastChar = nrChars-4 # .hdr / .img cutted
				if (index+3<=nrChars-4)	{
					for (i in (index+3):(nrChars-4)){
						mypostfix = paste(mypostfix,help2[i],sep="")
					}
				}
				startWith = as.numeric(help2[index])*100+as.numeric(help2[index+1])*10+as.numeric(help2[index+2])
				if (scanspS==1){
					#data <<- read.ANALYZE(dataFile) # hier koennte man noch 4D Analyze einbauen !!!
					dataFileNew = ""
					for (i in 1:(nrChars-4)){
						dataFileNew = paste(dataFileNew,help2[i],sep="")
					}
					data <<- read.ANALYZE(dataFileNew)	
				}
				else {
					data <<- read.ANALYZE(prefix=myprefix,postfix=mypostfix,numbered=TRUE,picstart=startWith,numbpic=scanspS)
				}
			}
			else {
				dataFileNew = ""
				for (i in 1:(nrChars-4)){
					dataFileNew = paste(dataFileNew,help2[i],sep="")
				}
				data <<- read.ANALYZE(dataFileNew)
			}
		}
		
		if (dataType=="NIFTI"){
			data <<- read.NIFTI(as.character(tclvalue(dataFileFirst)))
		}
		
				
		if (dataType=="unknown"){
			print("The data type is unknown !!")
			print("Please check your stated path or press 'help'.")
			tkconfigure(objFileLoad,bg="#FF0000")
		}
		else {
			print(paste("Data loaded  ...  data type: ",dataTypeGlobal,sep=""))
			if (data$dim[4]==scanspS){
				print("Check conformity of data and design  ...  Ok")
				tkconfigure(objFileLoad,bg=wiaslightblue)	
				conform <<- 1
			}	
			else {
				print("Check conformity of data and design ... Don't Conform !!")
				print("Please redefine the design or check the data.")
				tkconfigure(objFileLoad,bg="#FF0000")
				conform <<- 0
			}		
			if (conform==1){
				tclvalue(quantileTc) = round(100*quantile(extract.data(data),0.75))/100
				tkgrid(maskLabel1,padx=10,pady=10,sticky="ew")
				tkgrid(frameMask1)
				tkgrid(maskLabel2,maskEntry,maskButton1,padx=10,pady=10)
				tkgrid(frameMask2)
				tkgrid(maskButton2,helpButton3,padx=15,pady=10,sticky="ew")
				tkgrid(frameMask4)			
				tkgrid(helpLabel5,sticky="ew")
				tkgrid(frameMask3,sticky="ew")
			}			
		}		
	}

	helpFunction1 <- function() { helpFunction(1) }
	helpFunction2 <- function() { helpFunction(2) }
	helpFunction3 <- function() { helpFunction(3) }
	helpFunction4 <- function() { helpFunction(4) }
	helpFunction5 <- function() { helpFunction(5) }
	helpFunction6 <- function() { helpFunction(6) }
	
	helpFunction <- function(i){	
			onQuit <- function(){
				tkdestroy(ttHelp)
			}
			ttHelp = tktoplevel(bg=wiasblue)
			tkwm.title(ttHelp, "Help")
			helpFrame1 = tkframe(ttHelp,bg=wiasblue)
			helpLabel = tklabel(helpFrame1,text=helptextVec[i],font="Arial 13",bg=wiasblue)
			helpFrame2 = tkframe(ttHelp,bg=wiasblue)	
			helpB1 = tkbutton(helpFrame2,text="Quit",command=onQuit)
			tkgrid(helpLabel)
			tkgrid(helpB1,padx=10,pady=10)
			tkgrid(helpFrame1)	
			tkgrid(helpFrame2)	
	}
		

	
	# set variables
	helptextVec = c("a","b","c","d","e","f")
	wiasblue = "#AACCDB"
	wiaslightblue = "#BBDDEC"
	nrStep <<-0	
	mycontrastTc <- tclVar()
	mycontrast <- c()
	quantileTc <- tclVar()
	designFile <- tclVar("")
	dataFileFirst <- tclVar("")
	pathName <- tclVar("~/data/myproject/")	
	interscanIntTc <- tclVar()
	interscanInt <- -1	
	scanspSTc <- tclVar()
	scanspS <- -1	
	nrCondTc <- tclVar("1")
	nrCond <- -1
	namesTc <- c()
	names <- c()	
	durationsTc <- c()
	durations <- c()		
	onsetsTc <- c()	
	onsets <- c()			
	hMax <- tclVar("0")
	rbDesign <- tclVar("unused")
	dataType <- ""	
	pvalue <- list()
	spm <- list()
	data <- list()
	dataTypeGlobal = ""
	isFigure <<- FALSE
	index <<- -1
	indexhelp <<- -1
	conform = 0

	# base GUI window (base.aws)
	if(.Platform$OS.type == "windows") flush.console()
	base.aws <- tktoplevel(bg=wiasblue)	
    	tkwm.title(base.aws, "FMRI Analysis - AWS")

	# mainframes
	frameDesign1   <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameDesign2   <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameDesign3 <- tkframe(base.aws,relief="groove",borderwidth=2,bg=wiasblue)	
	frameData1     <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameData2     <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameData4    <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)	
	frameData3    <- tkframe(base.aws,relief="groove",borderwidth=2,bg=wiasblue)
	frameMask1 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameMask2 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameMask3 <- tkframe(base.aws,relief="groove",borderwidth=2,bg=wiasblue)
	frameMask4 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameContrast1 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameContrast2 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameContrast3    <- tkframe(base.aws,relief="groove",borderwidth=2,bg=wiasblue)
	frameResults <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameSmoothChoice1   <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameSmoothChoice2   <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)
	frameSmoothChoice3 <- tkframe(base.aws,relief="groove",borderwidth=2,bg=wiasblue)
	frameSmoothChoice4 <- tkframe(base.aws,relief="groove",borderwidth=0,bg=wiasblue)				
	intScLabel <- list()
	nrScLabel <- list()
	nrCoLabel <- list()	
	condLabel <- list()
	onsdurLabel <- list()	
	designChoiceLabel <- list()
	dataLabel1 <- list()
	dataLabel2 <- list()
	dataLabel3 <- list()
								
	# buttons, labels, entries 	
	objDesignL  <- tklabel(frameDesign1,text="Experimental Design",width=32,bg=wiasblue,font="Arial 13 bold")
	objDesignB1 <- tkbutton(frameDesign2,text="Load",width=15,command=loadDesignHelp,bg=wiaslightblue)
	objDesignB2 <- tkbutton(frameDesign2,text="Create",width=15,command=createDesignHelp,bg=wiaslightblue)
	helpButton1 <- tkbutton(frameDesign2,text="Help",width=15,command=helpFunction1,bg=wiaslightblue)	
	objFileL    <- tklabel(frameData1,text="Load Data",bg=wiasblue,font="Arial 13 bold")	
	objFileE1   <- tkentry(frameData2, textvariable = dataFileFirst, width = 40, bg = "#ffffff")
    	objFileB1   <- tkbutton(frameData2, text = "Select file", width = 15, command = selectFirstDataFile, bg = wiaslightblue, anchor = "c")	
	objFileLoad <- tkbutton(frameData4, text = "Load", width = 15, command = startAdjustHelp, bg = wiaslightblue)
	helpButton2 <- tkbutton(frameData4,text="Help",width=15,command=helpFunction2,bg=wiaslightblue)		
	maskLabel1  <- tklabel(frameMask1,text="Adjust Mask",bg=wiasblue,font="Arial 13 bold")	
	maskLabel2  <- tklabel(frameMask2,text="Threshold", width = 15, bg=wiasblue,font="Arial 12 bold")
	maskEntry   <- tkentry(frameMask2, textvariable = quantileTc, width = 15, bg = "#ffffff")
    	maskButton1 <- tkbutton(frameMask2, text = "View Mask", width = 20,command = viewMask, bg = wiaslightblue, anchor = "c")	
	maskButton2 <- tkbutton(frameMask4, text = "Ok", width = 15, command = startContrastHelp, bg = wiaslightblue, anchor = "c")	
	helpButton3 <- tkbutton(frameMask4,text="Help",width=15,command=helpFunction3,bg=wiaslightblue)	
	objcontrastL<- tklabel(frameContrast1,text="Define Contrast",bg=wiasblue,font="Arial 13 bold")
	objcontrastE<- tkentry(frameContrast2,textvariable = mycontrastTc, width = 20, bg = "#ffffff")
	objcontrastB<- tkbutton(frameContrast2,text="Ok",width=15,command=startEstimationHelp,bg=wiaslightblue)
	helpButton4 <- tkbutton(frameContrast2,text="Help",width=15,command=helpFunction4,bg=wiaslightblue)	
	objResultsL  <- tklabel(frameResults,text="View results",bg=wiasblue,font="Arial 12 bold")
	objResultsB1 <- tkbutton(frameResults,text="in slices",width=11,command=resultsSlices,bg=wiaslightblue)
	objResultsB2 <- tkbutton(frameResults,text="in 3D",width=11,command=results3D,bg=wiaslightblue)
	helpButton6 <- tkbutton(frameResults,text="Help",width=11,command=helpFunction6,bg=wiaslightblue)		
	objhMaxL <- tklabel(frameSmoothChoice2,text="Maximal bandwidth",bg=wiasblue,font="Arial 11 bold")
	objhMaxE <- tkentry(frameSmoothChoice2,textvariable=hMax,width=6,bg="#ffffff")
	objButtonSmooth <- tkbutton(frameSmoothChoice2,text="Start adaptive smoothing",width=20,command=startSmoothingHelp,bg=wiaslightblue)
	helpButton5 <- tkbutton(frameSmoothChoice4,text="Help",width=15,command=helpFunction5,bg=wiaslightblue)		
	helpLabel1 <- tklabel(frameData3,text="",bg=wiasblue,width=0,font="Arial 1")
	helpLabel2 <- tklabel(frameContrast3,text="",bg=wiasblue,width=0,font="Arial 1")
	helpLabel3 <- tklabel(frameDesign3,text="",bg=wiasblue,width=0,font="Arial 1")
	helpLabel4 <- tklabel(frameSmoothChoice3,text="",bg=wiasblue,width=0,font="Arial 1")
	helpLabel5 <- tklabel(frameMask3,text="",bg=wiasblue,width=0,font="Arial 1")	
	smoothChoiceL  <- tklabel(frameSmoothChoice1,text="Adaptive Smoothing",bg=wiasblue,font="Arial 13 bold")
	smoothChoiceB2 <- tkbutton(frameSmoothChoice4,text="Continue without adaptive smoothing",width=33,command=contWithoutHelp,bg=wiaslightblue)
				
	if (as.numeric(tclRequire("Img",warn=FALSE))!=0){
		image1 <- tclVar()
		tkimage.create("photo",image1,file=system.file("img/wias.jpeg",package="fmri"))
		imgAsLabel <- tklabel(base.aws,image=image1,bg="white")
		tkgrid(objDesignL,imgAsLabel,padx=20,pady=10)
	}
	else {
		tkgrid(objDesignL,pady = 10, padx = 10)
	}
	tkgrid(frameDesign1)	
	tkgrid(objDesignB1,objDesignB2,helpButton1,padx=15,pady=10)
	tkgrid(frameDesign2,sticky="ew")
 	tkgrid(helpLabel3)	
	tkgrid(frameDesign3,sticky="ew")

}
