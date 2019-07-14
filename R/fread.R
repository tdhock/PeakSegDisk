fread.first <- function
### Read the first line of a text file.
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  dt <- fread(file.name, nrows=1L, col.names=col.name.vec)
  dt
### Data table with one row.
}

fread.last <- function
### Read the last line of a text file.
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  wc.cmd <- paste("wc -l", file.name)
  wc.output <- system(wc.cmd, intern=TRUE)
  lines.chr <- sub(" .*", "", wc.output)
  lines.int <- as.integer(lines.chr)
  dt <- fread(file.name, skip=lines.int-1L, col.names=col.name.vec)
  dt
### Data table with one row.
}

