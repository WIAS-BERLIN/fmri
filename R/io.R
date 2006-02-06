read.ANALYZE <- function(prefix = "", numbered = FALSE, postfix = "", picstart = 1, numbpic = 1) {
  counter <- c(paste("00", 1:9, sep=""), paste("0", 10:99, sep=""),paste(100:999,sep=""));
  if (require(AnalyzeFMRI)) {
    if (numbered) {
      filename <- paste(prefix, counter[picstart], postfix, ".img", sep="")
    } else {
      filename <- paste(prefix, ".img", sep="")
    }
    
    if (length(system(paste("ls",filename),TRUE,TRUE)) != 0) {
      ttt <- f.read.analyze.volume(filename);
      dt <- dim(ttt)
      cat(".")
      header <- f.read.analyze.header(filename);

      if ((numbpic > 1) && !numbered) { 
        for (i in (picstart+1):(picstart+numbpic-1)) {
          filename <- paste(prefix, counter[i], postfix, ".img", sep="")
          a <- f.read.analyze.volume(filename)
          if (sum() != 0)
            cat("Error: wrong spatial dimension in picture",i)
          ttt <- c(ttt,a);
          dt[4] <- dt[4] + dim(a)[4]
          cat(".")
        }
      }

      cat("\n")
      dim(ttt) <- dt
      invisible(list(ttt=ttt,header=header))
    } else {
        warning(paste("Error: file",filename,"does not exist!"))
        list(ttt=NULL,header=NULL)
    }
  } else {
    warning("Error: library AnalyzeFMRI not found\n")
    list(ttt=NULL,header=NULL)
  }
}



write.ANALYZE <- function(ttt, name = "data", size="int", voxelsize = c(2,2,2)) {
  if (require(AnalyzeFMRI)) {
    f.write.analyze(ttt, file=name, size=size, voxelsize)
  } else {
    warning("Error: library AnalyzeFMRI not found\n")
    invisible(NULL)
  }
}



read.AFNI <- function(file) {
  conhead <- file(paste(file,".HEAD",sep=""),"r")
  header <- readLines(conhead)
  close(conhead)

  types <- NULL
  args <- NULL
  counts <- NULL
  values <- NULL
  
  for (i in 1:length(header)) {
    if (regexpr("^type *= *", header[i]) != -1) {
      tmptype <- strsplit(header[i]," *= *")[[1]][2]
      types <- c(types,tmptype)
      args <- c(args,strsplit(header[i+1]," *= *")[[1]][2])
      tmpcounts <- as.numeric(strsplit(header[i+2]," *= *")[[1]][2])
      counts <- c(counts,tmpcounts)
      i <- i+3
      tmpvalue <- ""
      while ((regexpr("^$", header[i]) == -1) && (i <= length(header))) {
        tmpvalue <- paste(tmpvalue,header[i])
        i <- i+1
      }
      tmpvalue <- sub("^ +","",tmpvalue)
      if ((tmptype == "integer-attribute") || (tmptype == "float-attribute")) {
        tmpvalue <- as.numeric(strsplit(tmpvalue," +")[[1]])
      }
      values <- c(values,list(value=tmpvalue))
    }        
  }

  names(values) <- args

  dx <- values$DATASET_DIMENSIONS[1]
  dy <- values$DATASET_DIMENSIONS[2]
  dz <- values$DATASET_DIMENSIONS[3]
  dt <- values$DATASET_RANK[2]
  scale <- values$BRICK_FLOAT_FACS
  size <- file.info(paste(file,".BRIK",sep=""))$size/(dx*dy*dz*dt)

  if (regexpr("MSB",values$BYTEORDER_STRING[1]) != -1) {
    endian <- "big"
  } else {
    endian <- "little"
  }

  if (as.integer(size) == size) {
    conbrik <- file(paste(file,".BRIK",sep=""),"rb")
    myttt<- readBin(conbrik, "int", n=dx*dy*dz*dt*size, size=size, signed=TRUE, endian=endian)
    close(conbrik)
    dim(myttt) <- c(dx,dy,dz,dt)
    for (k in 1:dt) {
      if (scale[k] != 0) {
        cat("scale",k,"with",scale[k],"\n")
        cat(range(myttt[,,,k]),"\n")
        myttt[,,,k] <- scale[k] * myttt[,,,k]
        cat(range(myttt[,,,k]),"\n")
      }
    }
    list(ttt=myttt,header=values)
  } else {
    warning("Error reading file: Could not detect size per voxel\n")
    list(ttt=NULL,header=values)    
  }
}




write.AFNI <- function(file, ttt, label, note="", origin=c(0,0,0), delta=c(4,4,4), idcode="WIAS_noid") {
  AFNIheaderpart <- function(type, name, value) {
    a <- "\n"
    a <- paste(a, "type = ", type, "\n", sep="")
    a <- paste(a, "name = ", name, "\n", sep="")
    if (regexpr("string",type) == 1) {
      value <- paste("'", value, "~", sep="")
      a <- paste(a, "count = ", nchar(value) - 1, "\n", sep ="")
      a <- paste(a, value, "\n", sep="")
    } else {
      a <- paste(a, "count = ", length(value), "\n", sep ="")
      j <- 0
      while (j<length(value)) {
        left <- length(value) - j
        if (left>4) left <- 5
        a <- paste(a, paste(value[(j+1):(j+left)],collapse="  "), "\n", sep="  ")
        j <- j+5
      }
    }
    a
  }
  
  conhead <- file(paste(file, ".HEAD", sep=""), "w")
  writeChar(AFNIheaderpart("string-attribute","HISTORY_NOTE",note),conhead,eos=NULL)
  writeChar(AFNIheaderpart("string-attribute","TYPESTRING","3DIM_HEAD_FUNC"),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_STRING",idcode),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_DATE",date()),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","SCENE_DATA",c(0,11,1,-999,-999,-999,-999,-999)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","ORIENT_SPECIFIC",c(0,3,4)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","ORIGIN",origin),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","DELTA",delta),conhead,eos=NULL)  
  minmax <- function(y) {r <- NULL;for (k in 1:dim(y)[4]) {r <- c(r,min(y[,,,k]),max(y[,,,k]))}; r}
  mm <- minmax(ttt)
  writeChar(AFNIheaderpart("float-attribute","BRICK_STATS",mm),conhead,eos=NULL)
  writeChar(AFNIheaderpart("integer-attribute","DATASET_RANK",c(3,dim(ttt)[4],0,0,0,0,0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","DATASET_DIMENSIONS",c(dim(ttt)[1:3],0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","BRICK_TYPES",rep(1,dim(ttt)[4])),conhead,eos=NULL)  

  scale <- rep(0,dim(ttt)[4])
  for (k in 1:dim(ttt)[4]) {
    scale[k] <- max(abs(mm[2*k-1]),abs(mm[2*k]))/32767
    ttt[,,,k] <- ttt[,,,k] / scale[k]
  }

  writeChar(AFNIheaderpart("float-attribute","BRICK_FLOAT_FACS",scale),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BRICK_LABS",paste(label,collapse="~")),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BYTEORDER_STRING","MSB_FIRST"),conhead,eos=NULL)  
  close(conhead)

  conbrik <- file(paste(file, ".BRIK", sep=""), "wb")
  dim(ttt) <- NULL
  writeBin(as.integer(ttt), conbrik,size=2, endian="big")
  close(conbrik)
}


