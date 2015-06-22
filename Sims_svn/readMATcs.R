###########################################################################/**
# \name{readMAT}
# \alias{readMAT}
#
# \title{Reads a MAT file structure from connection or file}
#
# \usage{
#   data <- readMAT(con, maxLength=NULL, verbose=FALSE)
# }
#
# \arguments{
#   \item{con}{Binary connection to which the MAT file structure should be
#     written to. A string is interpreted as filename, which then will be
#     opened (and closed afterwards).}
#   \item{maxLength}{The maximum number of bytes to be read from the input
#     stream, which should be equal to the length of the MAT file structure.
#     If \code{NULL}, data will be read until End Of File has been reached.
#     Default value is \code{NULL}.}
#   \item{verbose}{If \code{TRUE}, debug information is written to 
#     standard output, otherwise not. Default value is \code{FALSE}.}
# }
#
# \description{
#  Reads a MAT file structure from the input stream, either until End of File
#  is detected or until \code{maxLength} bytes has been read.
#  Using \code{maxLength} it is possible to read MAT file structure over
#  socket connections and other non-terminating input streams. In such cases
#  the \code{maxLength} has to be communicated before sending the MAT file
#  structure.
#
#  Currently only the MAT version 5 file format are supported.
# }
#
# \value{
#   Returns a named \code{list} structure containing all variables in the
#   MAT file structure.
# }
#
# \examples{
#  # See help on writeMAT() for readMAT() examples.
# }
#
# \author{Henrik Bengtsson, \url{http://www.braju.com/R/}}
#
# \seealso{
#   \code{\link{writeMAT}()}.
# }
#*/###########################################################################
readMAT <- function(con, maxLength=NULL, verbose=FALSE) {
  nbrOfBytesRead <- 0;
  detectedEndian <- "little";
  left <- NA;

  # Opens file if filename was given
  if (is.character(con)) {
    con <- file(con, "rb");
    on.exit(close(con));
  }

  # Assert that it is a binary connection
  if (summary(con)$text != "binary")
    stop("Can only write a MAT file structure to a *binary* connection.");


  willRead <- function(nbrOfBytes) {
  	if (is.null(maxLength))
  		return();
  	if (nbrOfBytesRead + nbrOfBytes <= maxLength)
  		return();
  	throw("Trying to read more bytes than expected from InputStream. Have read ", nbrOfBytesRead, " byte(s) and trying to read another ", nbrOfBytes, " byte(s), but expected ", maxLength, " byte(s).");
  } # willRead()
  
  
  hasRead <- function(nbrOfBytes) {
  	nbrOfBytesRead <<- nbrOfBytesRead + nbrOfBytes;
  	if (is.null(maxLength))
  		return(TRUE);
  	return(nbrOfBytesRead <= maxLength);
  } # hasRead()
  
  
  isDone <- function() {
  	if (is.null(maxLength))
  		return(FALSE);
  	return(nbrOfBytesRead >= maxLength);
  } # isDone()


  readBinMAT <- function(what, size, n, endian=detectedEndian) {
  	# Check maxLength to see if we are done.
  	if (isDone())
  		return(c());
  
  	willRead(size*n);
  	bfr <- readBin(con=con, what=what, size=size, n=n, endian=endian);
  	hasRead(length(bfr)*size);
  	bfr;
  } # readBinMAT()
  
  
  readCharMAT <- function(nchars) {
  	# Check maxLength to see if we are done.
  	if (isDone())
  		return(c());
  
  	willRead(nchars);
  	bfr <- readChar(con=con, nchars=nchars);
  	hasRead(nchars);
  	bfr;
  } # readCharMAT()


  readHeader <- function(this) {
  	# - - - - - - - - - - 
  	#  Text
  	# - - - - - - - - - - 
  	MOPT <- readBinMAT(what=integer(), size=1, n=4);
  	if (MOPT[1] %in% 0:4 && MOPT[2] == 0 && MOPT[3] %in% 0:5 && MOPT[4] %in% 0:2)
  		throw("MAT file format v4 is not supported.");
  	
  	description <- c(MOPT, readBinMAT(what=integer(), size=1, n=120));
  	description <- paste(intToChar(description), collapse="");
  
  	# - - - - - - - - - - 
  	#  Version
  	# - - - - - - - - - - 
  	# At this point we can not know which the endian is and we just have to
  	# make a guess and adjust later.  
  	version <- readBinMAT(what=integer(), size=2, n=1, endian="little");
  
  	# - - - - - - - - - - 
  	#  Endian Indicator
  	# - - - - - - - - - - 
  	endian <- readCharMAT(nchars=2);
  	if (endian == "MI")
  		detectedEndian <<- "big"
  	else if (endian == "IM")
  		detectedEndian <<- "little"
  	else {
  		warning(paste("Unknown endian: ", endian, ". Will assume Bigendian.", sep=""));
  		detectedEndian <<- "big";
  	}
  
  	if (detectedEndian == "big") {
  		 hi <- version %/% 256;
  		 low <- version %% 256;
  		 version <- 256*low + hi;
  	}
  
  	if (version == 256) {
  		version = "5";
  	} else {
  		warning(paste("Unknown MAT version tag: ", version, ". Will assume version 5.", sep=""));
  		version = as.character(version);
  	}
  
  	list(description=description, version=version);
  } # readHeader()


  readDataElement <- function(this) {
  	isSigned <- function(type) {
  		signed   <- c("mxINT8_CLASS", "mxINT16_CLASS", "mxINT32_CLASS");
  		signed   <- c(signed, "miINT8", "miINT16", "miINT32");
  		unsigned <- c("mxUINT8_CLASS", "mxUINT16_CLASS", "mxUINT32_CLASS");
  		unsigned <- c(unsigned, "miUINT8", "miUINT16", "miUINT32");
  		if (!is.element(type, c(signed, unsigned)))
  			return(NA);
  		is.element(type, signed);
  	} # isSigned()
  
  
  
  	readTag <- function(this) {
  		#      1    2    3    4    5    6    7    8
  		#   +----+----+----+----+----+----+----+----+
  		#   |    Data type      |  Number of Bytes  |  Tag
  		#   +---------------------------------------+
  		#   :                                       :
  		#
  		#   or...
  		#
  		#      1    2    3    4   ...
  		#   +----+----+----+----+----+----+----+----+
  		#   | Nbr of. | Data t. | ...               |  Tag
  		#   +---------------------------------------+
  		#   :                                       :
  		type <- readBinMAT(what=integer(), size=4, n=1);
  		# Did we read EOF?
  		if (length(type) == 0)
  			return(NULL);
  	
  		left <<- left - 4;
  		
  		knownTypes <- c("miINT8"=8, "miUINT8"=8, "miINT16"=16, "miUINT16"=16, "miINT32"=32, "miUINT32"=32, "miSINGLE"=NA, NA, "miDOUBLE"=64, NA, NA, "miINT64"=64, "miUINT64"=64, "miMATRIX"=NA);
  		knownWhats <- list("miINT8"=integer(), "miUINT8"=integer(), "miINT16"=integer(), "miUINT16"=integer(), "miINT32"=integer(), "miUINT32"=integer(), "miSINGLE"=NA, NA, "miDOUBLE"=double(), NA, NA, "miINT64"=integer(), "miUINT64"=integer(), "miMATRIX"=NA);
  	
  		compressed <- FALSE;
  		nbrOfBytes <- NULL;
  		if (!is.element(type, seq(knownTypes))) {
  			# Compressed tag?
  			nbrOfBytes <- type %/% 2^16;
  			type <- type %% 2^16;
  			# NOTE: Do not swap for different endians here. /HB 020827
  			if (!is.element(type, seq(knownTypes)))
  				throw("Unknown data type tag: ", type);
  			# Treat unsigned values too.
  			compressed <- TRUE;
  			padding <- 4 - ((nbrOfBytes-1) %% 4 + 1);
  		}
  	
  		type <- names(knownTypes)[type];
  		sizeOf <- as.integer(knownTypes[type]);
  		what <- knownWhats[[type]];
  	#  cat("type=", type, ", sizeOf=", sizeOf, ", what=", typeof(what), "\n", sep="");
  	
  		signed <- isSigned(type);
  		
  		if (!compressed) {
  			nbrOfBytes <- readBinMAT(what=integer(), size=4, n=1);
  			left <<- left - 4;
  			padding <- 8 - ((nbrOfBytes-1) %% 8 + 1);
  		}
  	
  		if (verbose)
  			cat("readTag(): ", type, ", ", nbrOfBytes, "\n", sep="");
  		list(type=type, signed=signed, sizeOf=sizeOf, what=what, nbrOfBytes=nbrOfBytes, padding=padding, compressed=compressed)
  	} # readTag()
  	
  	
  	readArrayFlags <- function(this) {
  		tag <- readTag(this);
  		
  		knownTypes <- c("mxCELL_CLASS"=NA, "mxSTRUCT_CLASS"=NA, "mxOBJECT_CLASS"=NA, "mxCHAR_CLASS"=8, "mxSPARSE_CLASS"=NA, "mxDOUBLE_CLASS"=NA, "mxSINGLE_CLASS"=NA, "mxINT8_CLASS"=8, "mxUINT8_CLASS"=8, "mxINT16_CLASS"=16, "mxUINT16_CLASS"=16, "mxINT32_CLASS"=32, "mxUINT32_CLASS"=32);
  	
  		arrayFlags <- readBinMAT(what=integer(), size=4, n=1);
  		left <<- left - 4;
  		
  		flags <- arrayFlags %/% 256;
  		flags <- as.logical(getBits(flags + 2^8)[-9]);
  		logical <- flags[2];
  		global  <- flags[3];
  		complex <- flags[4];
  		class <- arrayFlags %% 256;
  		if (!is.element(class, seq(knownTypes)))
  			throw("Unknown array type tag: ", class);
  		class <- names(knownTypes)[class];
  		classSize <- knownTypes[class];
  		
  		undefined <- readBinMAT(what=integer(), size=4, n=1);
  		left <<- left - 4;
  		
  		signed <- isSigned(tag$type);
  		
  		if (verbose)
  			cat("readArrayFlags(): ", class, "\n", sep="");
  		list(tag=tag, logical=logical, global=global, complex=complex, class=class, classSize=classSize, signed=signed);
  	} # readArrayFlags()
  	
  	
  	readDimensionsArray <- function(this) {
  		tag <- readTag(this);
  	
  		sizeOf <- tag$sizeOf %/% 8;
  		len <- tag$nbrOfBytes %/% sizeOf;
  		dim <- readBinMAT(what=integer(), size=sizeOf, n=len);
  		left <<- left - sizeOf*len;
  		
  		padding <- readBinMAT(what=integer(), size=1, n=tag$padding);
  		left <<- left - tag$padding;
  		
  		if (verbose)
  			cat("readDimensionsArray(): c(", paste(dim, collapse=","), ")\n", sep="");
  		list(tag=tag, dim=dim);
  	} # readDimensionsArray()
  
  
  	readName <- function(this) {
  		tag <- readTag(this);
  	
  		sizeOf <- tag$sizeOf %/% 8;
  		nchars <- tag$nbrOfBytes %/% sizeOf;
  		name <- readCharMAT(nchars=nchars);
  		left <<- left - nchars;
  	
  		padding <- readBinMAT(what=integer(), size=1, n=tag$padding);
  		left <<- left - tag$padding;
  		
  		if (verbose)
  			cat("readName(): '", name, "'\n", sep="");
  		list(tag=tag, name=name);
  	} # readName()
  	
  	
  	readFieldNameLength <- function(this) {
  		tag <- readTag(this);
  	
  		sizeOf <- tag$sizeOf %/% 8;
  		len <- tag$nbrOfBytes %/% sizeOf;
#  		cat("sizeOf=", sizeOf, "\n");
#  		cat("len=", len, "\n");
  		maxLength <- readBinMAT(what=integer(), size=sizeOf, n=len);
  		
  		left <<- left - len;
  	
  		padding <- readBinMAT(what=integer(), size=1, n=tag$padding);
  		left <<- left - tag$padding;
  	
  		if (verbose)
  			cat("readFieldNameLength(): ", maxLength, " field(s)\n", sep="");
  		list(tag=tag, maxLength=maxLength);
  	} # readFieldNameLength()
  	
  	
  	readFieldNames <- function(this, maxLength) {
  		tag <- readTag(this);
  	
  		names <- c();
  		nbrOfNames <- tag$nbrOfBytes %/% maxLength;
#  		cat("tag$nbrOfBytes=",tag$nbrOfBytes,"\n");
#  		cat("maxLength=",maxLength,"\n");
#  		cat("nbrOfNames=",nbrOfNames,"\n");
  		for (k in seq(nbrOfNames)) {
  	#    name <- readCharMAT(nchars=maxLength);
  			name <- readBinMAT(what=integer(), size=1, n=maxLength);
  			name <- intToChar(name);
  			name <- paste(name, collapse="");
  			left <<- left - maxLength;
  			names <- c(names, name);
  		}
  		
  		if (verbose)
  			cat("readFieldNames(): ", paste(paste("'", names, "'", sep=""), collapse=", "), "\n", sep="");
  		list(tag=tag, names=names);
  	} # readFieldNames()
  
  	readFields <- function(this, names) {
  		fields <- list();
  		for (k in seq(names)) {
  			field <- readDataElement(this);
  			fields <- c(fields, field);
  		}
  		names(fields) <- names;
  		
  		fields;
  	} # readFields()
  
  
  	readValues <- function(this) {
  		tag <- readTag(this);
  		sizeOf <- tag$sizeOf %/% 8;
  		len <- tag$nbrOfBytes %/% sizeOf;
  	
  		value <- readBinMAT(what=tag$what, size=sizeOf, n=len);
  		left <<- left - sizeOf*len;
  		
  		padding <- readBinMAT(what=integer(), size=1, n=tag$padding);
  		left <<- left - tag$padding;
  		
  		list(tag=tag, value=value);
  	} # readValues()
  
  
  
  	readMiMATRIX <- function(this) {
  		arrayFlags <- readArrayFlags(this);
  	#  print(arrayFlags);
  		dimensionsArray <- readDimensionsArray(this);
  	#  print(dimensionsArray);
  		arrayName <- readName(this);
  	#  print(arrayName);
  	
  		if (arrayFlags$class == "mxCELL_CLASS") {
  			nbrOfCells <- prod(dimensionsArray$dim);
  	#    cat("Reading ", nbrOfCells, " cells.\n", sep="");
  			matrix <- list();
  			for (k in seq(nbrOfCells)) {
  				tag <- readTag(this);
  				cell <- readMiMATRIX(this);
  				matrix <- c(matrix, cell);
  			}
  			matrix <- list(matrix);
  			names(matrix) <- arrayName$name;
  		} else if (arrayFlags$class == "mxSTRUCT_CLASS") {
  			maxLength <- readFieldNameLength(this);
  			names <- readFieldNames(this, maxLength=maxLength$maxLength);
  	#    print(names);
  			fields <- readFields(this, names=names$names);
  	#    str(fields);
  			matrix <- list(fields);
  			names(matrix) <- arrayName$name;
  		} else if (arrayFlags$class == "mxOBJECT_CLASS") {
  			className <- readName(this)$name;
  			maxLength <- readFieldNameLength(this);
  			names <- readFieldNames(this, maxLength=maxLength$maxLength);
  			fields <- readFields(this, names=names$names);
  			class(fields) <- className;
  			matrix <- list(fields);
  			names(matrix) <- arrayName$name;
  		} else if (arrayFlags$complex) {
  			pr <- readValues(this);
  			if (left > 0)
  				pi <- readValues(this);
  			matrix <- complex(real=pr$value, imaginary=pi$value);
  			dim(matrix) <- dimensionsArray$dim;
  			attr(matrix, "name") <- arrayName$name;
  		} else {
  			data <- readValues(this);
  			matrix <- data$value;
  	
  			if (arrayFlags$class == "mxDOUBLE_CLASS") {
  				matrix <- as.double(matrix);
  				dim(matrix) <- dimensionsArray$dim;
  			} else if (arrayFlags$class == "mxSINGLE_CLASS") {
  				matrix <- as.single(matrix);
  				dim(matrix) <- dimensionsArray$dim;
  			} else if (is.element(arrayFlags$class, c("mxINT8_CLASS", "mxUINT8_CLASS", "mxINT16_CLASS", "mxUINT16_CLASS", "mxINT32_CLASS", "mxUINT32_CLASS"))) {
  				matrix <- as.integer(matrix);
  				dim(matrix) <- dimensionsArray$dim;
  			} else if (arrayFlags$class == "mxCHAR_CLASS") {
  				matrix <- intToChar(matrix);
  				dim(matrix) <- dimensionsArray$dim;
  				matrix <- apply(matrix, MARGIN=1, FUN=paste, collapse="");
  				matrix <- as.matrix(matrix);
  			}
  	
  			matrix <- list(matrix);
  			names(matrix) <- arrayName$name;
  		}
  	
  		matrix;
  	} # readMiMATRIX()
  
  	#      1    2    3    4    5    6    7    8
  	#   +----+----+----+----+----+----+----+----+
  	#   |    Data type      |  Number of Bytes  |  Tag
  	#   +---------------------------------------+
  	#   |                                       |
  	#   |             Variable size             |  Data
  	#   |                                       |
  	#   +---------------------------------------+
  	tag <- readTag(this);
  	if (is.null(tag))
  		return(NULL);
  
  	left <<- tag$nbrOfBytes;
  	if (tag$type == "miMATRIX") {
  		data <- readMiMATRIX(this);
  	} else {
  		data <- readBinMAT(what=integer(), size=1, n=tag$nbrOfBytes);
  	}
  	data;
  } # readDataElement()


  nbrOfBytesRead <- 0;
  
  header <- readHeader(this);

  result <- list();
  data <- NA;
  while (!is.null(data)) {
    data <- readDataElement(this);
    result <- c(result, data);
  }
  result;
} # readMAT()


######################################################################
# HISTORY:
# 2002-09-03
# o readMAT() is now a stand-alone function.
# o Internal code need to be cleanup in the same fashion as in 
#   writeMAT.R.
# o Made readMAT() out of old MATInputStream.R.
# 2002-08-28
# o BUG FIX: The bug fix from yesterday where I thought the flag bits
#   should be readArrayFlags() in oposite order, was actually
#   incorrect. Excluded rev() again.
# 2002-08-27
# o TO DO: This class should be cleaned up in the same way as
#   MATOutputStream is.
# o BUG FIX: in readArrayFlags() the bits were read in the reverse
#   order, which resulted in incorrect values of complex, global, and
#   logical.
# o Updated readHeader() to read the endian information and possible
#   adjust the version tag if it was Bigendian. First, then the 
#   verification of the correct version number is done.
# o TEST: Went to the web and downloaded a few sample MAT files.
#   All v5 files loaded with any problems. The v4 files were reported
#   to be non-supported via an exception.
# o Added support for argument maxLength in read() and added
#   readMaxLength() for standardization (32-bit Littleendian integer)
#   and simplification.
# o Improved as.character().
# o Added some Rdoc comments.
# o Created.
######################################################################
