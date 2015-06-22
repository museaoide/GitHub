###########################################################################/**
# @RdocDefault readMat
#
# @title "Reads a MAT file structure from a connection or a file"
#
# \description{
#  Reads a MAT file structure from an input stream, either until End of File
#  is detected or until \code{maxLength} bytes has been read.
#  Using \code{maxLength} it is possible to read MAT file structure over
#  socket connections and other non-terminating input streams. In such cases
#  the \code{maxLength} has to be communicated before sending the actual
#  MAT file structure.
#
#  Both the MAT version 4 and MAT version 5 file formats are 
#  supported. The implementation is based on [1].
#
#  From Matlab v7, \emph{compressed} MAT version 5 files are used 
#  by default [3]. These are not supported. 
#  Use \code{save -V6} in Matlab to write a MAT file compatible with 
#  Matlab v6, that is, to write a non-compressed MAT version 5 file.
#  Note: Do not mix up version numbers for the Matlab software and
#  the Matlab file formats.
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{Binary @connection to which the MAT file structure should be
#     written to. A string is interpreted as filename, which then will be
#     opened (and closed afterwards).}
#   \item{maxLength}{The maximum number of bytes to be read from the input
#     stream, which should be equal to the length of the MAT file structure.
#     If        f \code{NULL}, data will be read until End Of File has been reached.}
#   \item{fixNames}{If @TRUE, names of Matlab variables and fields are 
#     renamed such that they are valid variables names in R.}
#   \item{verbose}{If @TRUE, debug information is written to standard output,
#     otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure containing all variables in the
#   MAT file structure.
# }
#
# \details{
#   For the MAT v5 format, \emph{cell} structures are read into
#   \R as a @list structure.
#
#   Sparse matrices are converted into plain matrices, which means that
#   some matrices will be too large to be allocated.
# }
#
# @examples "readMat.Rex"
#
# \author{
#   Henrik Bengtsson, Mathematical Statistics, Lund University.
#   The internal MAT v4 reader was written by 
#   Andy Jacobson at Program in Atmospheric and Oceanic Sciences, 
#   Princeton University. 
# }
#
# \seealso{
#   @see "writeMat".
# }
#
# \references{
#   [1] The MathWorks Inc., \emph{Matlab - MAT-File Format, version 5}, June 1999.\cr
#   [2] The MathWorks Inc., \emph{Matlab - Application Program Interface Guide, version 5}, 1998.\cr
#   [3] The MathWorks Inc., \emph{Matlab - MAT-File Format, version 7}, October 2004, \url{http://www.mathworks.com/access/helpdesk/help/pdf_doc/matlab/matfile_format.pdf}\cr
# }
#
# @keyword file
# @keyword IO
#*/###########################################################################
#setMethodS3("readMat", "default", function(con, maxLength=NULL, fixNames=TRUE, verbose=FALSE, ...) {
readMat <- function(con, maxLength=NULL, fixNames=TRUE, verbose=FALSE, ...) {
                                        #===========================================================================
                                        # General functions to read both MAT v4 and MAT v5 files.              BEGIN
                                        #===========================================================================
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # willRead(), hasHead() and isDone() operators keep count on the number of
                                        # bytes actually read and compares it with 'maxLength'.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  nbrOfBytesRead <- 0;

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # readBinMat() need to know what endian the numerics in the stream are 
                                        # written with. From the beginning we assume Little Endian, but that might
                                        # be updated when we have read the MAT-file header.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  detectedEndian <- "little"; 

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # ASCII is the 8-bit ASCII table with ASCII characters from 0-255.
                                        # 
                                        # Extracted from the R.oo package.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ASCII <- c(
             "\000","\001","\002","\003","\004","\005","\006","\007", # 000-007
             "\010","\011","\012","\013","\014","\015","\016","\017", # 010-017
             "\020","\021","\022","\023","\024","\025","\026","\027", # 020-027
             "\030","\031","\032","\033","\034","\035","\036","\037", # 030-037
             "\040","\041","\042","\043","\044","\045","\046","\047", # 040-047
             "\050","\051","\052","\053","\054","\055","\056","\057", # 050-057
             "\060","\061","\062","\063","\064","\065","\066","\067", # 060-067
             "\070","\071","\072","\073","\074","\075","\076","\077", # 070-077
             "\100","\101","\102","\103","\104","\105","\106","\107", # 100-107
             "\110","\111","\112","\113","\114","\115","\116","\117", # 110-117
             "\120","\121","\122","\123","\124","\125","\126","\127", # 120-127
             "\130","\131","\132","\133","\134","\135","\136","\137", # 130-137
             "\140","\141","\142","\143","\144","\145","\146","\147", # 140-147
             "\150","\151","\152","\153","\154","\155","\156","\157", # 150-157
             "\160","\161","\162","\163","\164","\165","\166","\167", # 160-167
             "\170","\171","\172","\173","\174","\175","\176","\177", # 170-177
             "\200","\201","\202","\203","\204","\205","\206","\207", # 200-207
             "\210","\211","\212","\213","\214","\215","\216","\217", # 210-217
             "\220","\221","\222","\223","\224","\225","\226","\227", # 220-227
             "\230","\231","\232","\233","\234","\235","\236","\237", # 230-237
             "\240","\241","\242","\243","\244","\245","\246","\247", # 240-247
             "\250","\251","\252","\253","\254","\255","\256","\257", # 250-257
             "\260","\261","\262","\263","\264","\265","\266","\267", # 260-267
             "\270","\271","\272","\273","\274","\275","\276","\277", # 270-277
             "\300","\301","\302","\303","\304","\305","\306","\307", # 300-307
             "\310","\311","\312","\313","\314","\315","\316","\317", # 310-317
             "\320","\321","\322","\323","\324","\325","\326","\327", # 320-327
             "\330","\331","\332","\333","\334","\335","\336","\337", # 330-337
             "\340","\341","\342","\343","\344","\345","\346","\347", # 340-347
             "\350","\351","\352","\353","\354","\355","\356","\357", # 350-357
             "\360","\361","\362","\363","\364","\365","\366","\367", # 360-367
             "\370","\371","\372","\373","\374","\375","\376","\377" # 370-377
             );

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to convert a vector of integers into a vector of ASCII chars.
                                        # 
                                        # Extracted from the R.oo package.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  intToChar <- function(i) {
    ASCII[i %% 256 + 1];
  } 

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to assert that it is possible to read a certain number of bytes.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  willRead <- function(nbrOfBytes) {
    if (is.null(maxLength))
      return();
    if (nbrOfBytesRead + nbrOfBytes <= maxLength)
      return();
    stop(paste("Trying to read more bytes than expected from connection. Have read ", nbrOfBytesRead, " byte(s) and trying to read another ", nbrOfBytes, " byte(s), but expected ", maxLength, " byte(s).", sep=""));
  }                                     # willRead()
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to tell now many bytes we actually have read.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  hasRead <- function(nbrOfBytes) {
    nbrOfBytesRead <<- nbrOfBytesRead + nbrOfBytes;
    if (is.null(maxLength))
      return(TRUE);
    return(nbrOfBytesRead <= maxLength);
  }                                     # hasRead()
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to check is there are more bytes to read.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  isDone <- function() {
    if (is.null(maxLength))
      return(FALSE);
    return(nbrOfBytesRead >= maxLength);
  }                                     # isDone()



                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read 'n' binary values of a data type of type 'what' 
                                        # and size 'size', cf. readBin(). 
                                        # This function will also keep track of the actual number of bytes read.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readBinMat <- function(con, what, size, n, signed=TRUE, endian=detectedEndian) {
                                        # Check maxLength to see if we are done.
    if (isDone())
      return(c());
    if (is.na(signed))
      signed <- TRUE;  
    willRead(size*n);
    bfr <- readBin(con=con, what=what, size=size, n=n, signed=signed, endian=endian);
                                        #    print(c(size=size, n=n, signed=signed, endian=endian, bfr=bfr));
    hasRead(length(bfr)*size);
    bfr;
  }                                     # readBinMat()
    
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read 'nchars' characters from the input connection.
                                        # This function will also keep track of the actual number of bytes read.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readCharMat <- function(con, nchars) {
                                        # Check maxLength to see if we are done.
    if (isDone())
      return(c());
  
    willRead(nchars);
    bfr <- readChar(con=con, nchars=nchars);
    hasRead(nchars);
    bfr;
  }                                     # readCharMat()


                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to make a variable name into a safe R variable name.
                                        # For instance, underscores ('_') are replaced by periods ('.').
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  asSafeRName <- function(name) {
    if (fixNames) {
      name <- gsub("_", ".", name);
    }
    name;
  }

  debugIndent <- 0;
  debug <- function(..., sep="") {
    cat(paste(rep(" ", length.out=debugIndent), collapse=""));
    cat(..., sep=sep);
    cat("\n");
  }

  debugPrint <- function(...) {
    print(...);
  }

  debugStr <- function(...) {
    str(...);
  }

  debugEnter <- function(..., indent=+1) {
    debug(..., "...");
    debugIndent <<- debugIndent + indent;
  }
  
  debugExit <- function(..., indent=-1) {
    debugIndent <<- debugIndent + indent;
    debug(..., "...done\n");
  }
  
                                        #===========================================================================
                                        # General functions to read both MAT v4 and MAT v5 files.                END
                                        #===========================================================================


                                        #===========================================================================
                                        # MAT v4 specific                                                      BEGIN
                                        #===========================================================================

                                        # "Programming Note When creating a MAT-file, you must write data in the
                                        #  first four bytes of this header. MATLAB uses these bytes to determine
                                        #  if a MAT-file uses a Version 5 format or a Version 4 format. If any of
                                        #  these bytes contain a zero, MATLAB will incorrectly assume the file is
                                        #  a Version 4 MAT-file."
  isMat4 <- function(MOPT) {
    any(MOPT == 0);
  }

  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to convert four signed or unsigned integers in big or little
                                        # endian order into a (MOPT) vector c(M,O,P,T) of unsigned integers.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  getMOPT <- function(fourBytes) {
    if (length(fourBytes) != 4)
      stop(paste("Argument 'fourBytes' must a vector of 4 bytes:", length(fourBytes)));
    
                                        # Make sure the four bytes are non-signed integers
    fourBytes <- as.integer(fourBytes);
    neg <- (fourBytes < 0);
    if (any(neg))
      fourBytes[neg] <- fourBytes[neg] + 256;

    base <- 256^(0:3);
    MOPT <- c(NA,NA,NA,NA);
    for (endian in c("little", "big")) {
      mopt <- sum(base*fourBytes);
      for (kk in 4:1) {
        MOPT[kk] <- mopt %% 10;
        mopt <- mopt %/% 10;
      }

      isMOPT <- (MOPT[1] %in% 0:4 && MOPT[2] == 0 && MOPT[3] %in% 0:5 && MOPT[4] %in% 0:2);
      if (isMOPT)
        break;

      base <- rev(base);
    }                                   # for (endian ...)

    if (!isMOPT)
      stop("File format error: Not a valid MAT v4. The first four bytes (MOPT) were: ", paste(MOPT, collapse=", "));
    
    if (verbose)
      cat("Read MOPT bytes: ", moptToString(MOPT), "\n", sep="");
    
    MOPT;
  }                                     # getMOPT()


 
  readMat4 <- function(con, maxLength=NULL, firstFourBytes=NULL) {
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read a MAT v4 Matrix Header Format
                                        # 
                                        # Fix length: 20 bytes
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    readMat4Header <- function(con, firstFourBytes=NULL) {
      header <- list();

                                        # "A MAT-file may contain one or more matrices. The matrices are written
                                        #  sequentially on disk, with the bytes forming a continuous stream. Each matrix
                                        #  starts with a fixed-length 20-byte header that contains information describing
                                        #  certain attributes of the Matrix. The 20-byte header consists of five long
                                        #  (4-byte) integers."


                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'type' field, a.k.a. MOPT
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(firstFourBytes)) {
        firstFourBytes <- readBinMat(con, what=integer(), size=1, n=4);
      }

                                        # If no bytes are read, we have reached the End Of Stream.
      if (length(firstFourBytes) == 0)
        return(NULL);

                                        # Assert that it really is a MAT v4 file we are reading and get MOPT bytes
      MOPT <- getMOPT(firstFourBytes);

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[1] "indicates the numeric format of binary numbers on the machine 
                                        #          that wrote the file.
                                        #          0 IEEE Little Endian (PC, 386, 486, DEC Risc)
                                        #          1 IEEE Big Endian (Macintosh, SPARC, Apollo,SGI, HP 9000/300,
                                        #            other Motorola)
                                        #          2 VAX D-float  [don't know how to read these]
                                        #          3 VAX G-float  [don't know how to read these]
                                        #          4 Cray         [don't know how to read these]"
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (MOPT[1] == 0) {
  	detectedEndian <<- "little";
      } else if (MOPT[1] == 1) {
  	detectedEndian <<- "big";
      } else if (MOPT[1] %in% 2:4) {
  	stop("Looks like a MAT v4 file, but the storage format of numerics (VAX D-float, VAX G-float or Cray) is not supported. Currently only IEEE numeric formats in big or little endian are supported.");
      } else {
  	stop(paste("Unknown first byte in MOPT header (not in [0,4]): ", paste(MOPT, collapse=", ")));
      }

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[2] "is always 0 (zero) and is reserved for future use."
                                        #
                                        # I've only seen a non-default value for ocode used once, by a
                                        # matfile library external to the MathWorks.  I believe it stands
                                        # for "order" code...whether a matrix is written in row-major or
                                        # column-major format.  Its value here will be ignored. /Andy November 2003
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      header$ocode <- MOPT[2];
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[3] "indicates which format the data is stored in according to the 
                                        #          following table:
                                        #          0 double-precision (64-bit) floating point numbers
                                        #          1 single-precision (32-bit) floating point numbers
                                        #          2 32-bit signed integers
                                        #          3 16-bit signed integers
                                        #          4 16-bit unsigned integers
                                        #          5 8-bit unsigned integers
                                        #          The precision used by the save command depends on the size and
                                        #          type of each matrix. Matrices with any noninteger entries and 
                                        #          matrices with 10,000 or fewer elements are saved in floating 
                                        #          point formats requiring 8 bytes per real element. Matrices 
                                        #          with all integer entries and more than 10,000 elements are 
                                        #          saved in the following formats, requiring fewer bytes per element."
                                        #
                                        # precision defines the number of type of data written and thus the number of 
                                        # bytes per datum.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (MOPT[3] == 0) {
                                        # "64-bit double";
        header$what <- double();
        header$size <- 8;
        header$signed <- NA;
      } else if (MOPT[3] == 1) {
                                        # "32-bit single";
        header$what <- double();
        header$size <- 4;
        header$signed <- NA;
      } else if (MOPT[3] == 2) {
                                        # "32-bit signed integer";
        header$what <- integer();
        header$size <- 4;
        header$signed <- TRUE;          # Ignored by readBin() because 32-bit ints are always signed!
      } else if (MOPT[3] == 3) {
                                        # "16-bit signed integer";
        header$what <- integer();
        header$size <- 2;
        header$signed <- TRUE;
      } else if (MOPT[3] == 4) {
                                        # "16-bit unsigned integer";
        header$what <- integer();
        header$size <- 2;
        header$signed <- FALSE;
      } else if (MOPT[3] == 5) {
                                        # "8-bit unsigned integer";
        header$what <- integer();
        header$size <- 1;
        header$signed <- FALSE;
      } else {
  	stop(paste("Unknown third byte in MOPT header (not in [0,5]): ", paste(MOPT, collapse=", ")));
      }
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[4]  "indicates the matrix type according to the following table:
                                        #          0 Numeric (Full) matrix
                                        #          1 Text matrix
                                        #          2 Sparse matrix
                                        #          Note that the elements of a text matrix are stored as floating
                                        #          point numbers between 0 and 255 representing ASCII-encoded 
                                        #          characters."
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      header$matrixType <- "numeric";
      if (MOPT[4] == 0) {
        header$matrixType <- "numeric";
      } else if (MOPT[4] == 1) {
        header$matrixType <- "text";
      } else if (MOPT[4] == 2) {
        header$matrixType <- "sparse";
      } else {
                                        #  	stop(paste("Unknown fourth byte in MOPT header (not in [0,2]): ", paste(MOPT, collapse=", ")));
      }

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'mrows' and 'ncols' fields
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "The row dimension contains an integer with the number of rows in the matrix."
      header$mrows  <- readBinMat(con, what=integer(), size=4, n=1)

                                        # "The column dimension contains an integer with the number of columns in the matrix."
      header$ncols  <- readBinMat(con, what=integer(), size=4, n=1)

      if (verbose)
        cat("Matrix dimension: ", header$mrows, "x", header$ncols, "\n", sep="");

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'imagf' fields
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "The imaginary flag is an integer whose value is either 0 or 1. If 1, 
                                        #  then the matrix has an imaginary part. If 0, there is only real data."
      header$imagf  <- readBinMat(con, what=integer(), size=4, n=1)

      if (verbose)
        cat("Matrix contains imaginary values: ", as.logical(header$imagf), "\n", sep="");

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'namelen' fields
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "The name length contains an integer with 1 plus the length of the matrix name."
      header$namlen <- readBinMat(con, what=integer(), size=4, n=1)

      if (verbose)
        cat("Matrix name length: ", header$namlen-1, "\n", sep="");

      header;
    }
    
    
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read a MAT v4 Matrix Data Format
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    readMat4Data <- function(con, header) {
                                        # "Immediately following the fixed length header is the data whose length
                                        #  is dependent on the variables in the fixed length header:"

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'name' field
                                        #
                                        # "The matrix name consists of 'namlen' ASCII bytes, the last one of which
                                        #  must be a null character (\0)."
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      name <- readCharMat(con, header$namlen);

      if (verbose)
        cat("Matrix name: '", name, "'\n", sep="");
      
      name <- asSafeRName(name);

      if (verbose)
        cat("Matrix safe name: '", name, "'\n", sep="");

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # The 'real' field
                                        #
                                        # "Real part of the matrix consists of mrows * ncols numbers in the format
                                        #  specified by the MOPT[3] element of the type flag. The data is stored
                                        #  column-wise such that the second column follows the first column, etc."
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n <- header$mrows * header$ncols;
      if (header$matrixType == "text") {
  	##data <- readCharMat(con, nchars=n);
        ##data <- strsplit(data, split="");
        data <- readBinMat(con,what=header$what,size=header$size,signed=header$signed, n=n)
        data <- as.integer(data)
        data <- intToChar(data)
                                        # Make into a matrix
        dim(data) <- c(header$mrows, header$ncols);
        data <- apply(data,MARGIN=1,FUN=function(s){ paste(s,sep="", collapse="")})
      } else if (header$matrixType %in% c("numeric", "sparse")) {
  	real <- readBinMat(con, what=header$what, size=header$size, signed=header$signed, n=n);
  	if (header$imagf != 0) {
  	  imag <- readBinMat(con, what=header$what, size=header$size, signed=header$signed, n=n);
  	  data <- complex(real=real, imag=imag);
  	} else {
  	  data <- real;
  	  rm(real);
  	}
        
                                        # Make into a matrix
        dim(data) <- c(header$mrows, header$ncols);
        
        if (header$matrixType == "sparse") {
                                        # From help sparse in Matlab:
                                        # "S = SPARSE(i,j,s,m,n,nzmax) uses the rows of [i,j,s] to generate an
                                        #  m-by-n sparse matrix with space allocated for nzmax nonzeros.  The
                                        #  two integer index vectors, i and j, and the real or complex entries
                                        #  vector, s, all have the same length, nnz, which is the number of
                                        #  nonzeros in the resulting sparse matrix S .  Any elements of s
                                        #  which have duplicate values of i and j are added together."
          i <- data[,1];
          j <- data[,2];
          s <- data[,3];
          rm(data);
          
                                        # When save a sparse matrix, Matlab is making sure that one can infer
                                        # the size of the m-by-n sparse matrix for the index matrix [i,j]. If
                                        # there are no non-zero elements in the last row or last column, Matlab
                                        # saves a zero elements in such case.
          n <- max(i);
          m <- max(j);
          
                                        # Instead of applying row-by-row, we calculate the position of each sparse
                                        # element in an hardcoded fashion.
          pos <- (j-1)*n + i;
          rm(i,j);                      # Not needed anymore

          pos <- pos[-length(pos)];
          s <- s[-length(s)];
          
          data <- matrix(0, nrow=n, ncol=m);
          data[pos] <- s;

          rm(pos, s);                   # Not needed anymore
        }
      } else {
        stop(paste("MAT v4 file format error: Unknown 'type' in header: ", header$matrixType, sep=""));
      }

      if (verbose) {
        cat("Matrix elements:\n");
        str(data);
      }
      
      data <- list(data);
      names(data) <- name;

      data;
    }
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "Main program"
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # Since readMat4() is wrapped inside the readMat() function, we can assume
                                        # that 'con' really is a connection.

    result <- list();
  
    repeat {
      header <- readMat4Header(con, firstFourBytes=firstFourBytes);
      if (is.null(header))
  	break;

      data <- readMat4Data(con, header);
      result <- append(result, data);
      rm(data);

      firstFourBytes <- NULL;
    }                                   # repeat

    header <- list(version="4", endian=detectedEndian);
    attr(result, "header") <- header;
    result;
  }                                     # readMat4()

  
                                        # Debug function to generate more informative error messages.
  moptToString <- function(MOPT) {
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[1] "indicates the numeric format of binary numbers on the machine 
                                        #          that wrote the file.
                                        #          0 IEEE Little Endian (PC, 386, 486, DEC Risc)
                                        #          1 IEEE Big Endian (Macintosh, SPARC, Apollo,SGI, HP 9000/300,
                                        #            other Motorola)
                                        #          2 VAX D-float  [don't know how to read these]
                                        #          3 VAX G-float  [don't know how to read these]
                                        #          4 Cray         [don't know how to read these]"
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (MOPT[1] == 0)
      mStr <- "IEEE Little Endian (PC, 386, 486, DEC Risc)"
    else if (MOPT[1] == 1)
      mStr <- "IEEE Big Endian (Macintosh, SPARC, Apollo,SGI, HP 9000/300, other Motorola)"
    else if (MOPT[1] == 2)
      mStr <- "VAX D-float"
    else if (MOPT[1] == 3)
      mStr <- "VAX G-float"
    else if (MOPT[1] == 4)
      mStr <- "Cray"
    else
      mStr <- sprintf("<Unknown value of MOPT[1]. Not in range [0,4]: %d.>", as.integer(MOPT[1]));

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[2] "is always 0 (zero) and is reserved for future use."
                                        #
                                        # I've only seen a non-default value for ocode used once, by a
                                        # matfile library external to the MathWorks.  I believe it stands
                                        # for "order" code...whether a matrix is written in row-major or
                                        # column-major format.  Its value here will be ignored. /Andy November 2003
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (MOPT[2] == 0)
      oStr <- "Reserved for future use"
    else
      oStr <- sprintf("<Unknown value of MOPT[2]. Should be 0: %d.>", as.integer(MOPT[2]));

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[3] "indicates which format the data is stored in according to the 
                                        #          following table:
                                        #          0 double-precision (64-bit) floating point numbers
                                        #          1 single-precision (32-bit) floating point numbers
                                        #          2 32-bit signed integers
                                        #          3 16-bit signed integers
                                        #          4 16-bit unsigned integers
                                        #          5 8-bit unsigned integers
                                        #          The precision used by the save command depends on the size and
                                        #          type of each matrix. Matrices with any noninteger entries and 
                                        #          matrices with 10,000 or fewer elements are saved in floating 
                                        #          point formats requiring 8 bytes per real element. Matrices 
                                        #          with all integer entries and more than 10,000 elements are 
                                        #          saved in the following formats, requiring fewer bytes per element."
                                        #
                                        # precision defines the number of type of data written and thus the number of 
                                        # bytes per datum.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (MOPT[3] == 0)
      pStr <- "64-bit double"
    else if (MOPT[3] == 1)
      pStr <- "32-bit single"
    else if (MOPT[3] == 2)
      pStr <- "32-bit signed integer"
    else if (MOPT[3] == 3)
      pStr <- "16-bit signed integer"
    else if (MOPT[3] == 4)
      pStr <- "16-bit unsigned integer"
    else if (MOPT[3] == 5)
      pStr <- "8-bit unsigned integer"
    else
      pStr <- sprintf("<Unknown value of MOPT[3]. Not in range [0,5]: %d.>", as.integer(MOPT[3]));

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # MOPT[4]  "indicates the matrix type according to the following table:
                                        #          0 Numeric (Full) matrix
                                        #          1 Text matrix
                                        #          2 Sparse matrix
                                        #          Note that the elements of a text matrix are stored as floating
                                        #          point numbers between 0 and 255 representing ASCII-encoded 
                                        #          characters."
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (MOPT[4] == 0)
      tStr <- "Numeric (Full) matrix"
    else if (MOPT[4] == 1)
      tStr <- "Text matrix"
    else if (MOPT[4] == 2)
      tStr <- "Sparse matrix"
    else
      tStr <- sprintf("<Unknown value of MOPT[4]. Not in range [0,2]: %d.>", as.integer(MOPT[4]));


    moptStr <- paste("MOPT[1]: ", mStr, ". MOPT[2]: ", oStr, ". MOPT[3]: ", pStr, ". MOPT[4]: ", tStr, ".", sep="");
    moptStr;
  }                                     # moptToString()
  
                                        #===========================================================================
                                        # MAT v4 specific                                                        END
                                        #===========================================================================

                                        #===========================================================================
                                        # MAT v5 specific                                                      BEGIN
                                        #===========================================================================
  readMat5 <- function(con, maxLength=NULL, firstFourBytes=NULL) {
                                        # Used to test if there a matrix read contains an imaginary part too.
    left <- NA;

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read the MAT-file header, which contains information of what
                                        # version of the MAT file we are reading, if it used little or big endian
                                        # etc.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    readMat5Header <- function(this, firstFourBytes=NULL) {
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "MATLAB uses the first four bytes to determine if a MAT-file uses a
                                        #  Version 5 format or a Version 4 format. If any of these bytes
                                        #  contain a zero, MATLAB will assume the file is a Version 4 MAT-file."
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(firstFourBytes))
  	firstFourBytes <- readBinMat(con, what=integer(), size=1, n=4);
  
      MOPT <- firstFourBytes;

      if (MOPT[1] %in% 0:4 && MOPT[2] == 0 && MOPT[3] %in% 0:5 && MOPT[4] %in% 0:2) {
  	stop("Detected MAT file format v4. Do not use readMat5() explicitly, but use readMat().");
      }
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        #  Text [124 bytes] (we already have read four of them)
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      description <- c(MOPT, readBinMat(con, what=integer(), size=1, n=120));
      description <- paste(intToChar(description), collapse="");
                                        #      cat("Description: '", description, "'\n", sep="");
  
                                        # - - - - - - - - - - 
                                        #  Version
                                        # - - - - - - - - - - 
                                        # At this point we can not know which the endian is and we just have to
                                        # make a guess and adjust later. 
      version <- readBinMat(con, what=integer(), size=2, n=1, endian="little");

                                        # - - - - - - - - - - 
                                        #  Endian Indicator
                                        # - - - - - - - - - - 
      endian <- readCharMat(con, nchars=2);
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
  
      if (version == 256) {             # version == 0x0100
  	version = "5";
      } else {
  	warning(paste("Unknown MAT version tag: ", version, ". Will assume version 5.", sep=""));
  	version = as.character(version);
      }
    
      list(description=description, version=version, endian=detectedEndian);
    }                                   # readMat5Header()
  
  
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                        # Function to read a MAT v5 Data Element
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    readMat5DataElement <- function(this) {
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      isSigned <- function(type) {
  	signed   <- c("mxINT8_CLASS", "mxINT16_CLASS", "mxINT32_CLASS");
  	signed   <- c(signed, "miINT8", "miINT16", "miINT32");
  	unsigned <- c("mxUINT8_CLASS", "mxUINT16_CLASS", "mxUINT32_CLASS");
  	unsigned <- c(unsigned, "miUINT8", "miUINT16", "miUINT32");
  	if (!is.element(type, c(signed, unsigned)))
  	  return(NA);
  	is.element(type, signed);
      }                                 # isSigned()
       
    
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
                                        #
                                        # "Each data element begins with an 8-byte tag followed immediately
                                        #  by the data in the element."
                                        #
                                        # From [1, page 6]:
                                        #
                                        #      1    2    3    4    5    6    7    8
                                        #   +----+----+----+----+----+----+----+----+
                                        #   |    Data type      |  Number of Bytes  |  Tag
                                        #   +---------------------------------------+
                                        #   :                                       :
                                        #
                                        # but also [1, page 9]:
                                        #
                                        #      1    2    3    4   ...
                                        #   +----+----+----+----+----+----+----+----+
                                        #   | Nbr of. | Data t. | ...               |  Tag
                                        #   +---------------------------------------+
                                        #   :                                       :
                                        #
      readTag <- function(this) {
        if (verbose)
          debugEnter("Reading Tag");
        
  	type <- readBinMat(con, what=integer(), size=4, n=1);
                                        # Did we read EOF?
  	if (length(type) == 0)
  	  return(NULL);
      
  	left <<- left - 4;

        knownTypes <- c("miMATRIX"=0, "miINT8"=8, "miUINT8"=8, "miINT16"=16, "miUINT16"=16, "miINT32"=32, "miUINT32"=32, "miSINGLE"=NA, NA, "miDOUBLE"=64, NA, NA, "miINT64"=64, "miUINT64"=64, "miMATRIX"=NA);
  	knownWhats <- list("miMATRIX"=0, "miINT8"=integer(), "miUINT8"=integer(), "miINT16"=integer(), "miUINT16"=integer(), "miINT32"=integer(), "miUINT32"=integer(), "miSINGLE"=NA, NA, "miDOUBLE"=double(), NA, NA, "miINT64"=integer(), "miUINT64"=integer(), "miMATRIX"=NA);
      
  	nbrOfBytes <- NULL;

                                        # From [1, page 9]:
                                        # "Programming Note - When reading a MAT-file, you can tell if you 
                                        #  are processing a compressed data element by comparing the value
                                        #  of the first two bytes of the tag with the value zero (0). If
                                        #  these two bytes are not zero, the tag uses the compressed format."
        tmp <- type;
        bytes <- rep(NA, length=4);
        for (kk in 1:4) {
          bytes[kk] <- (tmp %% 256);
          tmp <- tmp %/% 256;
        }
        rm(tmp);
        compressed <- any(bytes[3:4] != 0);

        if (verbose)
          debug("Compressed tag: ", compressed);

                                        #  	  if (type+1 < 1 || type+1 > length(knownTypes)) {
        if (compressed) {
                                        #           stop()
                                        # NOTE: Do not swap for different endians here. /HB 020827
  	  nbrOfBytes <- type %/% 2^16;
  	  type <- type %% 2^16;
          if (detectedEndian == "big") {
            tmp <- type;
                                        #            type <- nbrOfBytes;
                                        #            nbrOfBytes <- tmp;
          }
  	  if (type+1 < 1 || type+1 > length(knownTypes))
  	    stop(paste("Unknown data type. Not in range [1,", length(knownTypes), "]: ", type, sep=""));
          
                                        # Treat unsigned values too.
  	  padding <- 4 - ((nbrOfBytes-1) %% 4 + 1);
  	} else {
                                        #    print(c(size=size, n=n, signed=signed, endian=endian, bfr=bfr));
  	  nbrOfBytes <- readBinMat(con, what=integer(), size=4, n=1);
  	  left <<- left - 4;
  	  padding <- 8 - ((nbrOfBytes-1) %% 8 + 1);
        }
      
  	type <- names(knownTypes)[type+1];
  	sizeOf <- as.integer(knownTypes[type]);
  	what <- knownWhats[[type]];
                                        #        cat("type=", type, ", sizeOf=", sizeOf, ", what=", typeof(what), "\n", sep="");
      
  	signed <- isSigned(type);
      
  	tag <- list(type=type, signed=signed, sizeOf=sizeOf, what=what, nbrOfBytes=nbrOfBytes, padding=padding, compressed=compressed);
        
        if (verbose) {
          debugPrint(unlist(tag));
          debugExit("Reading Tag");
        }
        
        tag;
      }                                 # readTag()
      
      
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
                                        #
                                        # Subelement     Data Type  Number of Bytes
                                        # ---------------------------------------------------------------------
                                        # Array Flags    miUINT32   2*sizeOf(miUINT32) (8 bytes)
      
      readArrayFlags <- function(this) {
        if (verbose)
          debugEnter("Reading Array Flags");
        
  	getBits <- function(i) {
  	  ready <- FALSE;
  	  bits <- c();
  	  while (!ready) {
  	    bit <- i %% 2;
  	    bits <- c(bits, bit);
  	    i <- i %/% 2;
  	    ready <- (i==0);
  	  }
  	  bits;
  	}                               # getBits()
  
  	tag <- readTag(this);
  
  	knownTypes <- c("mxCELL_CLASS"=NA, "mxSTRUCT_CLASS"=NA, "mxOBJECT_CLASS"=NA, "mxCHAR_CLASS"=8, "mxSPARSE_CLASS"=NA, "mxDOUBLE_CLASS"=NA, "mxSINGLE_CLASS"=NA, "mxINT8_CLASS"=8, "mxUINT8_CLASS"=8, "mxINT16_CLASS"=16, "mxUINT16_CLASS"=16, "mxINT32_CLASS"=32, "mxUINT32_CLASS"=32);

                                        # Read the first miUINT32 integer
  	arrayFlags <- readBinMat(con, what=integer(), size=4, n=1);
  	left <<- left - 4;

                                        # Byte 4 - Class
                                        # "Class. This field contains a value that identifies the MATLAB
                                        # array type (class) represented by the data element."
                                        #
  	class <- arrayFlags %% 256;
  	if (class < 1 || class > length(knownTypes)) { 
  	  stop(paste("Unknown array type (class). Not in [1,",
                     length(knownTypes), "]: ", class, sep=""));
        }
 	class <- names(knownTypes)[class];
  	classSize <- knownTypes[class];

       	arrayFlags <- arrayFlags %/% 256;

                                        # Byte 3 - Flags
                                        # "Flags. This field contains three, single-bit flags that indicate
                                        #  whether the numeric data is complex, global, or logical. If the
                                        #  complex bit is set, the data element includes an imaginary part
                                        #  (pi). If the global bit is set, MATLAB loads the data element as
                                        #  a global variable in the base workspace. If the logical bit is
                                        #  set, it indicates the array is used for logical indexing."
  	flags <- arrayFlags %% 256;
  	flags <- as.logical(getBits(flags + 2^8)[-9]);
  	logical <- flags[2];
  	global  <- flags[3];
  	complex <- flags[4];

                                        # Bytes 1 & 2 - The two hi-bytes are "undefined".


                                        # Used for Sparse Arrays, otherwise undefined
                                        # Read the second miUINT32 integer
  	nzmax <- readBinMat(con, what=integer(), size=4, n=1);
  	left <<- left - 4;
  	
  	signed <- isSigned(tag$type);
  	
  	flags <- list(tag=tag, logical=logical, global=global, complex=complex, class=class, classSize=classSize, signed=signed, nzmax=nzmax);

        if (verbose) {
          debugPrint(unlist(flags[-1]));
          debugExit("Reading Array Flags");
        }
        
        flags;
      }                                 # readArrayFlags()
      
      
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readDimensionsArray <- function(this) {
        if (verbose)
          debugEnter("Reading Dimensions Array");
        
  	tag <- readTag(this);
      
  	sizeOf <- tag$sizeOf %/% 8;
  	len <- tag$nbrOfBytes %/% sizeOf;
        if (verbose)
          debug("Reading ", len, " integers each of size ", sizeOf, " bytes.");
  	dim <- readBinMat(con, what=integer(), size=sizeOf, n=len);
  	left <<- left - sizeOf*len;
  	
        if (verbose)
          debug("Reading ", tag$padding, " padding bytes.");
  	padding <- readBinMat(con, what=integer(), size=1, n=tag$padding);
  	left <<- left - tag$padding;
  	
  	dimArray <- list(tag=tag, dim=dim);
        
        if (verbose) {
          debugPrint(list(dim=dim));
          debugExit("Reading Dimensions Array");
        }
        
        dimArray;
      }                                 # readDimensionsArray()
    
    
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readName <- function(this) {
        if (verbose)
          debugEnter("Reading Array Name");
        
  	tag <- readTag(this);
      
  	sizeOf <- tag$sizeOf %/% 8;
  	nchars <- tag$nbrOfBytes %/% sizeOf;
        if (verbose)
          debug("Reading ", nchars, " characters.");
  	name <- readCharMat(con, nchars=nchars);
        name <- asSafeRName(name);
  	left <<- left - nchars;
      
        if (verbose)
          debug("Reading ", tag$padding, " padding bytes.");
  	padding <- readBinMat(con, what=integer(), size=1, n=tag$padding);
  	left <<- left - tag$padding;
  	
        if (verbose) {
  	  debug("Name: '", name, "'");
          debugExit("Reading Array Name");
        }
        
  	list(tag=tag, name=name);
      }                                 # readName()
      
      
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readFieldNameLength <- function(this) {
        if (verbose)
          debugEnter("Reading Field Name Length");
        
  	tag <- readTag(this);
      
  	sizeOf <- tag$sizeOf %/% 8;
  	len <- tag$nbrOfBytes %/% sizeOf;
                                        #      cat("sizeOf=", sizeOf, "\n");
                                        #      cat("len=", len, "\n");
  	maxLength <- readBinMat(con, what=integer(), size=sizeOf, n=len);
  	
  	left <<- left - len;
      
  	padding <- readBinMat(con, what=integer(), size=1, n=tag$padding);
  	left <<- left - tag$padding;
      
  	if (verbose)
  	  debug("Field name length+1: ", maxLength);

        if (verbose)
          debugExit("Reading Field Name Length");
        
  	list(tag=tag, maxLength=maxLength);
      }                                 # readFieldNameLength()
      
      
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readFieldNames <- function(this, maxLength) {
        if (verbose)
          debugEnter("Reading Field Names");
        
  	tag <- readTag(this);
      
  	names <- c();
  	nbrOfNames <- tag$nbrOfBytes %/% maxLength;
                                        #      cat("tag$nbrOfBytes=",tag$nbrOfBytes,"\n");
                                        #      cat("maxLength=",maxLength,"\n");
                                        #      cat("nbrOfNames=",nbrOfNames,"\n");
  	for (k in seq(length=nbrOfNames)) {
                                        #    name <- readCharMat(con, nchars=maxLength);
  	  name <- readBinMat(con, what=integer(), size=1, n=maxLength);
  	  name <- intToChar(name);
  	  name <- paste(name, collapse="");
          name <- asSafeRName(name);
  	  left <<- left - maxLength;
  	  names <- c(names, name);
  	}

        if (verbose)
          debug("Reading ", tag$padding, " padding bytes.");
  	padding <- readBinMat(con, what=integer(), size=1, n=tag$padding);
  	left <<- left - tag$padding;

  	if (verbose)
  	  debug("Field names: ", paste(paste("'", names, "'", sep=""), collapse=", "));
        
        if (verbose)
          debugExit("Reading Field Names");
        
  	list(tag=tag, names=names);
      }                                 # readFieldNames()
    

                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
                                        # From [1, page 26]:
                                        # "Fields Subelement - This subelement contains the value stored in a
                                        #  field. These values are MATLAB arrays, represented using the
                                        #  miMATRIX format specific to the array type: numeric array, sparse
                                        #  array, cell, object or other structure. See the appropriate section
                                        #  of this document for details about the MAT-file format of each of
                                        #  these array type. MATLAB reads and writes these fields in
                                        #  column-major order."
      readFields <- function(this, names) {
        if (verbose)
          debugEnter("Reading Fields");
        
  	fields <- list();
  	for (k in seq(names)) {
          if (verbose)
            debugEnter("Reading field: ", names[k]);
  	  field <- readMat5DataElement(this);
  	  fields <- c(fields, field);
          if (verbose)
            debugExit("Reading field: ", names[k]);
  	}
  	names(fields) <- names;
  	
        if (verbose)
          debugExit("Reading Fields");
        
  	fields;
      }                                 # readFields()
    
    
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readValues <- function(this) {
        if (verbose)
          debugEnter("Reading Values");
        
  	tag <- readTag(this);
  	sizeOf <- tag$sizeOf %/% 8;
  	len <- tag$nbrOfBytes %/% sizeOf;

        if (verbose)
          debug("Reading ", len, " values each of ", sizeOf, " bytes. In total ", tag$nbrOfBytes, " bytes.");
        
  	value <- readBinMat(con, what=tag$what, size=sizeOf, n=len, signed=tag$signed);
        if (verbose)
          debugStr(value);
        
  	left <<- left - sizeOf*len;
  	
        if (verbose)
          debug("Reading ", tag$padding, " padding bytes.");
        
  	padding <- readBinMat(con, what=integer(), size=1, n=tag$padding);
        
  	left <<- left - tag$padding;
  	
        if (verbose)
          debugExit("Reading Values");
        
  	list(tag=tag, value=value);
      }                                 # readValues()
    
    
    
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      readMiMATRIX <- function(this) {
        if (verbose)
          debugEnter("Reading miMATRIX");
        
  	arrayFlags <- readArrayFlags(this);
  	dimensionsArray <- readDimensionsArray(this);
  	arrayName <- readName(this);

                                        #str(arrayName)
      
  	if (arrayFlags$class == "mxCELL_CLASS") {
  	  nbrOfCells <- prod(dimensionsArray$dim);
          if (verbose)
            cat("Reading ", nbrOfCells, " cells.\n", sep="");
  	  matrix <- list();
  	  for (kk in seq(length=nbrOfCells)) {
  	    tag <- readTag(this);
  	    cell <- readMiMATRIX(this);
  	    matrix <- c(matrix, cell);
  	  }
  	  matrix <- list(matrix);
  	  names(matrix) <- arrayName$name;
  	} else if (arrayFlags$class == "mxSTRUCT_CLASS") {
  	  nbrOfCells <- prod(dimensionsArray$dim);
          if (verbose)
            cat("Reading ", nbrOfCells, " cells in structure.\n", sep="");
  	  maxLength <- readFieldNameLength(this);
  	  names <- readFieldNames(this, maxLength=maxLength$maxLength);
          if (verbose)
            cat("Field names: ", paste(names$names, collapse=", "), "\n", sep="");
          nbrOfFields <- length(names$names);
  	  matrix <- list();
  	  for (kk in seq(length=nbrOfCells)) {
                                        #            cat("Cell: ", kk, "...\n", sep="");
    	    fields <- readFields(this, names=names$names);
                                        #            str(fields);
            matrix <- c(matrix, fields);
                                        #            cat("Cell: ", kk, "...done\n", sep="");
  	  }
          names(matrix) <- NULL;

                                        # Set the dimension of the structure
          dim <- c(nbrOfFields, dimensionsArray$dim);
          if (prod(dim) > 0) {
            matrix <- structure(matrix, dim=dim);
            dimnames <- rep(list(NULL), length(dim(matrix)));
            dimnames[[1]] <- names$names;
            dimnames(matrix) <- dimnames;
          }

                                        # Finally, put the structure in a named list.
          matrix <- list(matrix);
          names(matrix) <- arrayName$name;

          if (verbose) {
            cat("Read a 'struct':\n");
            str(matrix);
          }
                                        #old#  	  maxLength <- readFieldNameLength(this);
                                        #old#  	  names <- readFieldNames(this, maxLength=maxLength$maxLength);
                                        #old#         print(names);
                                        #old#  	  fields <- readFields(this, names=names$names);
                                        #old#  	  names(fields) <- arrayName$name;
                                        #old#          str(fields);
                                        #old#  	  matrix <- list(fields);
                                        #old#  	  names(matrix) <- arrayName$name;
  	} else if (arrayFlags$class == "mxOBJECT_CLASS") {
  	  className <- readName(this)$name;
  	  maxLength <- readFieldNameLength(this);
  	  names <- readFieldNames(this, maxLength=maxLength$maxLength);
  	  fields <- readFields(this, names=names$names);
  	  class(fields) <- className;
  	  matrix <- list(fields);
  	  names(matrix) <- arrayName$name;
          ## Lines here are fix from Bengtsson, June 05
          } else if (arrayFlags$complex) {
            pr <- readValues(this);
            if (left > 0)
              pi <- readValues(this);

            matrix <- complex(real=pr$value, imaginary=pi$value);
            ## NOT NEEDED: attr(matrix, "name") <- arrayName$name;

                                        # Set dimension of complex matrix
            dim(matrix) <- dimensionsArray$dim; ## NEW

                                        # Put into a named list
            matrix <- list(matrix);     ## NEW
            names(matrix) <- arrayName$name; ## NEW
          } else if (arrayFlags$class == "mxSPARSE_CLASS") { 
            ##
          ##} else if (arrayFlags$complex) {
          ##  pr <- readValues(this);
          ##  if (left > 0)
          ##    pi <- readValues(this);
          ##  matrix <- complex(real=pr$value, imaginary=pi$value);
          ##  attr(matrix, "name") <- arrayName$name;
          ##} else if (arrayFlags$class == "mxSPARSE_CLASS") {
                                        # Dimensions of the sparse matrix
            nrow <- dimensionsArray$dim[1];
            ncol <- dimensionsArray$dim[2];
          
                                        # Create expanded matrix...
            matrix <- matrix(0, nrow=nrow, ncol=ncol);
            attr(matrix, "name") <- arrayName$name;
          
                                        # From [2, page5-6]
                                        # "Sparse Matrices
                                        #  Sparse matrices have a different storage convention in MATLAB. The
                                        #  parameters pr and pi are still arrays of double-precision numbers,
                                        #  but there are three additional parameters, nzmax, ir, and jc:
                                        #  * nzmax - is an integer that contains the length of ir, pr, and,
                                        #    if it exists, pi. It is the maximum possible number of nonzero
                                        #    elements in the sparse matrix.
                                        #  * ir - points to an integer array of length nzmax containing the
                                        #    row indices of the corresponding elements in pr and pi.
                                        #  * jc - points to an integer array of length N+1 that contains
                                        #    column index information. For j, in the range 0 <= j <= N-1,
                                        #    jc[j] is the index in ir and pr (and pi if it exists) of the
                                        #    first nonzero entry in the jth column and jc[j+1] - 1 index of
                                        #    the last nonzero entry. As a result, jc[N] is also equal to nnz,
                                        #    the number of nonzero entries in the matrix. If nnz is less
                                        #    than nzmax, then more nonzero entries can be inserted in the
                                        #    array without allocating additional storage."

                                        # From mxGetIr in [2]:
                                        # "The nzmax field holds an integer value that signifies the number
                                        #  of elements in the ir, pr, and, if it exists, the pi arrays. The
                                        #  value of nzmax is always greater than or equal to the number of
                                        #  nonzero elements in a sparse mxArray. In addition, the value of
                                        #  nzmax is always less than or equal to the number of rows times
                                        #  the number of columns."
            nzmax <- arrayFlags$nzmax;
            if (nzmax > 0) {
                                        # Read the row indices for non-zero values (index start at zero!)
                                        #
                                        # From mxGetIr in [2]:
                                        # "Each value in an ir array indicates a row (offset by 1) at which
                                        #  a nonzero element can be found. (The jc array is an index that
                                        #  indirectly specifies a column where nonzero elements can be found.)"
                                        #
                                        # and from mxSetIr in [2]:
                                        # "The ir array must be in column-major order. That means that the
                                        #  ir array must define the row positions in column 1 (if any) first,
                                        #  then the row positions in column 2 (if any) second, and so on through
                                        #  column N. Within each column, row position 1 must appear prior to
                                        #  row position 2, and so on."
              ir <- readValues(this)$value;
              if (length(ir) != nzmax) {
                stop(paste("MAT v5 file format error: The length of row index vector 'ir' (sparse arrays) is not equal to 'nzmax': ", length(ir), ", ", nzmax, "."));
              }
            
                                        # Note that the indices for MAT v5 sparse arrays start at 0 (not 1).
              ir <- ir + 1;
              if (any(ir < 1 | ir > nrow)) {
                stop(paste("MAT v5 file format error: Some elements in row vector 'ir' (sparse arrays) are out of range [1,", nrow, "].", sep=""));
              }

                                        #  "* jc - points to an integer array of length N+1 that contains..."
              jc <- readValues(this)$value;
              if (length(jc) != ncol+1) {
                stop(paste("MAT v5 file format error: Length of column vector 'jc' (sparse arrays) is not ", ncol, "+1 as expected: ", length(jc)));
              }
                                        # Add one to all indices except last one
              jc <- jc + 1;
              jc[length(jc)] <- jc[length(jc)] + 1;

                                        # Read real part
              pr <- readValues(this)$value;
              if (length(pr) != nzmax) {
                stop(paste("MAT v5 file format error: The length of vector 'pr' (sparse arrays) is not equal to 'nzmax': ", length(ir), ", ", nzmax, "."));
              }
  
              if (verbose) {
                debugStr(ir);
                debugStr(jc);
                debugStr(pr);
              }
            
                                        # "This subelement contains the imaginary data in the array, if one
                                        #  or more of the numeric values in the MATLAB array is a complex
                                        #  number (if the complex bit is set in Array Flags)." [1, p20]
              if (arrayFlags$complex) {
                                        # Read imaginary part
                pi <- readValues(this)$value;
                if (length(pi) != nzmax) {
                  stop(paste("MAT v5 file format error: The length of vector 'pi' (sparse arrays) is not equal to 'nzmax': ", length(ir), ", ", nzmax, "."));
                }
                if (verbose)
                  debugStr(pi);
                pr <- complex(real=pr, imaginary=pi);
                rm(pi);                 # Not needed anymore!
              }
           
                                        # Now, for each column insert the non-zero elements
                                        #
                                        #  "* jc - points to an integer array of length N+1 that contains
                                        #    column index information. For j, in the range 0 <= j <= N-1,
                                        #    jc[j] is the index in ir and pr (and pi if it exists) of the
                                        #    first nonzero entry in the jth column and jc[j+1] - 1 index of
                                        #    the last nonzero entry. As a result, jc[N] is also equal to nnz,
                                        #    the number of nonzero entries in the matrix. If nnz is less
                                        #    than nzmax, then more nonzero entries can be inserted in the
                                        #    array without allocating additional storage."
                                        #
                                        #    Note: This is *not* how MAT v4 works.
              for (col in seq(length=length(jc)-1)) {
                first <- jc[col];
                last  <- jc[col+1]-1;
                idx <- seq(from=first, to=last);
                value <- pr[idx];
                row <- ir[idx];
                ok <- is.finite(row);
                row <- row[ok];
                value <- value[ok];
                matrix[row,col] <- value;
              }
              rm(ir,jc,first,last,idx,value,row); # Not needed anymore
            
              matrix <- list(matrix);
              names(matrix) <- arrayName$name;
            }
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
            } else {
              stop(paste("Unknown or unsupported class id in array flags: ", arrayFlags$class, sep=""));
            }
      
            matrix <- list(matrix);
            names(matrix) <- arrayName$name;
          }
      
        if (verbose)
          debugExit("Reading miMATRIX");
        
  	matrix;
      }                                 # readMiMATRIX()
    

                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
                                        # General structure of a MAT v5 Data Element:
                                        # 
                                        #      1    2    3    4    5    6    7    8
                                        #   +----+----+----+----+----+----+----+----+
                                        #   |    Data type      |  Number of Bytes  |  Tag
                                        #   +---------------------------------------+
                                        #   |                                       |
                                        #   |             Variable size             |  Data
                                        #   |                                       |
                                        #   +---------------------------------------+
                                        #
                                        # -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
      tag <- readTag(this);
      if (is.null(tag))
  	return(NULL);
    
      if (tag$nbrOfBytes == 0)
  	return(list(NULL));

      left <<- tag$nbrOfBytes;
      if (tag$type == "miMATRIX") {
  	data <- readMiMATRIX(this);
      } else {
  	data <- readBinMat(con, what=integer(), size=1, n=tag$nbrOfBytes, signed=tag$signed);
      }

      data;
    }                                   # readMat5DataElement()

      
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "Main program"
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # Since readMat5() is wrapped inside the readMat() function, we can
                                        # assume that 'con' really is a connection.

    detectedEndian <<- "little";

    header <- readMat5Header(this, firstFourBytes=firstFourBytes);

    if (verbose) {
      debug("Read MAT v5 header:");
      debugPrint(header);
      debug("Endian: ", detectedEndian);
    }

    result <- list();
    repeat {
      data <- readMat5DataElement(this);
      if (is.null(data))
        break;
      result <- append(result, data);
    }

    attr(result, "header") <- header;
    result;
  }                                     # readMat5()
                                        #===========================================================================
                                        # MAT v5 specific                                                        END
                                        #===========================================================================

  if (inherits(con, "connection")) {
    if (!isOpen(con)) {
      open(con, open="rb");
      on.exit(close(con));
    }
  } else {
                                        # For all other types of values of 'con' make it into a character string.
                                        # This will for instance also make it possible to use object of class
                                        # File in the R.io package to be used.
    con <- as.character(con);

                                        # Now, assume that 'con' is a filename specifying a file to be opened.
    con <- file(con, open="rb");
    on.exit(close(con));
  }

                                        # Assert that it is a binary connection that we are reading from
  if (summary(con)$text != "binary")
    stop("Can only read a MAT file structure from a *binary* connection.");

                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                        # "MATLAB uses the first four bytes to determine if a MAT-file uses a
                                        #  Version 5 format or a Version 4 format. If any of these bytes
                                        #  contain a zero, MATLAB will assume the file is a Version 4 MAT-file."
                                        #
                                        # Details:
                                        # For MAT v5 the first 124 bytes is free text followed by 2 bytes 
                                        # specifying the version and 2 bytes specifying the endian, but for 
                                        # MAT v4 the first four bytes represents the 'type'. 
                                        #
                                        # Thus, we read the first four bytes and test if it can be a MAT v4 file.
                                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfBytesRead <- 0;
  firstFourBytes <- readBinMat(con, what=integer(), size=1, n=4);
  if (is.null(firstFourBytes))
    stop("MAT file format error: Nothing to read. Empty input stream.");

  if (isMat4(firstFourBytes)) {
    readMat4(con, firstFourBytes=firstFourBytes, maxLength=maxLength);
  } else {
    readMat5(con, firstFourBytes=firstFourBytes, maxLength=maxLength);
  }
}
