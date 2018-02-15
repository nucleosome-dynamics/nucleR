#####  Tools for reading paired ended data
####################################################################
readSingleFile <- function (fname, dir=getwd(), type='Bowtie', paired=FALSE,
                            ...)
{
    if (file.exists(file.path(dir,fname))) {
        cat('Reading ',fname,'...\n',sep='')
        pattern <- paste("^",fname,"$",sep='')
        ans <- readAligned(dirPath=dir,pattern=pattern,type=type,...)
        if (paired) {
            if (type == 'Bowtie') {
                cat(' Processing paired reads...')
                xplus <- ans[seq(1, length(ans), by=2)]
                xminus <- ans[seq(2, length(ans), by=2)]
                chrnames <- levels(ans@chromosome)
                ans <- vector('list', length(chrnames))
                names(ans) <- chrnames
                for (j in 1:length(chrnames)) {
                    sel <- which(chromosome(xplus) == chrnames[j])
                    start <- position(xplus)[sel],
                    end <- position(xminus)[sel] + width(xminus)[sel]
                    ans[[j]] <- IRanges(start=start, end=end)
                }
                ans <- RangedData(IRangesList(ans))
            }
        } else {
            if (length(grep('^chr', ans@chromosome)) == 0) {
                ans@chromosome <- factor(paste('chr', ans@chromosome, sep=''))
            }
            ans <- RangedData(ranges = IRanges(start = position(ans),
                                               width = width(ans)),
                              strand = strand(ans),
                              space  = ans@chromosome)
        }
        gc()
  } else {
      cat(paste('File',file.path(dir,fname),'does not exist\n'))
      ans <- NA
  }
  return(ans)
}

readAlignedBatch <- function (fnames, dir=getwd(), type='Bowtie', paired=FALSE,
                              mc.cores=1, ...)
{
  if (mc.cores != 1) {
      require(multicore)
  }
  if (paired & type!='Bowtie') {
      warning('Cannot deal with paired information for data aligned with ',
              type,
              '. Putting all reads into a single object.')
  }
  fnameslist <- vector("list", length(fnames))
  for (i in 1:length(fnameslist)) {
      fnameslist[[i]] <- fnames[i]
  }
  if (mc.cores == 1) {
      ans <- lapply(fnameslist,
                    readSingleFile,
                    dir    = dir,
                    type   = type,
                    paired = paired,
                    ...)
  } else {
      ans <- mclapply(fnameslist,
                      readSingleFile,
                      dir    = dir,
                      type   = type,
                      paired = paired,
                      ...,
                      mc.preschedule = FALSE,
                      mc.cores       = mc.cores)
  }
  ans <- RangedDataList(ans)
  names(ans) <- sub("\\..+$", '', fnames)
  return(ans)
}
