# Minimal ABIF reader for Sanger .ab1 files.
# Reads only tags required by this app: FWO_1, PBAS1/2, PCON1/2, PLOC1/2, DATA9-12.

raw_u8 <- function(x) as.numeric(as.integer(x))

be_u16 <- function(x) {
  b <- raw_u8(x)
  if (length(b) != 2L) stop("Expected 2 bytes")
  b[1] * 256 + b[2]
}

be_i16_vec <- function(x) {
  b <- raw_u8(x)
  if (length(b) %% 2L) stop("Odd byte count for int16")
  m <- matrix(b, ncol=2, byrow=TRUE)
  v <- m[,1] * 256 + m[,2]
  v[v >= 32768] <- v[v >= 32768] - 65536
  as.numeric(v)
}

be_u16_vec <- function(x) {
  b <- raw_u8(x)
  if (length(b) %% 2L) stop("Odd byte count for uint16")
  m <- matrix(b, ncol=2, byrow=TRUE)
  as.numeric(m[,1] * 256 + m[,2])
}

be_u32 <- function(x) {
  b <- raw_u8(x)
  if (length(b) != 4L) stop("Expected 4 bytes")
  b[1] * 256^3 + b[2] * 256^2 + b[3] * 256 + b[4]
}

be_i32_vec <- function(x) {
  b <- raw_u8(x)
  if (length(b) %% 4L) stop("Bad byte count for int32")
  m <- matrix(b, ncol=4, byrow=TRUE)
  v <- m[,1]*256^3 + m[,2]*256^2 + m[,3]*256 + m[,4]
  v[v >= 2^31] <- v[v >= 2^31] - 2^32
  as.numeric(v)
}

raw_slice0 <- function(raw, offset0, n) {
  if (n <= 0) return(raw())
  i <- as.integer(offset0) + 1L
  j <- i + as.integer(n) - 1L
  if (i < 1L || j > length(raw)) stop("ABIF offset outside file")
  raw[i:j]
}

parse_dir_entry <- function(x) {
  if (length(x) != 28L) stop("ABIF directory entry must be 28 bytes")
  list(
    tag=rawToChar(x[1:4]),
    number=be_u32(x[5:8]),
    element_type=be_u16(x[9:10]),
    element_size=be_u16(x[11:12]),
    n_elements=be_u32(x[13:16]),
    data_size=be_u32(x[17:20]),
    data_offset=be_u32(x[21:24]),
    inline=x[21:24],
    handle=be_u32(x[25:28])
  )
}

entry_bytes <- function(file_raw, e) {
  n <- as.integer(e$data_size)
  if (n <= 4L) e$inline[seq_len(n)] else raw_slice0(file_raw, e$data_offset, n)
}

entry_value <- function(file_raw, e) {
  b <- entry_bytes(file_raw, e)
  typ <- as.integer(e$element_type)
  # ABIF types: 1 byte, 2 char, 3 word, 4 short, 5 long, 7 float, 8 double.
  if (typ == 1L) return(raw_u8(b))
  if (typ == 2L) {
    txt <- rawToChar(b, multiple=FALSE)
    return(sub("\\x00+$", "", txt))
  }
  if (typ == 3L) return(be_u16_vec(b))
  if (typ == 4L) return(be_i16_vec(b))
  if (typ == 5L) return(be_i32_vec(b))
  # For tags this app needs, other types are not expected. Preserve raw bytes.
  b
}

read_abif <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add=TRUE)
  n <- file.info(path)$size
  r <- readBin(con, what="raw", n=n)
  if (length(r) < 34L || rawToChar(r[1:4]) != "ABIF") stop("Not an ABIF/AB1 file: ", basename(path))
  version <- be_u16(r[5:6])
  root <- parse_dir_entry(r[7:34])
  n_dir <- as.integer(root$n_elements)
  dir_off <- as.integer(root$data_offset)
  if (n_dir < 1L || n_dir > 100000L) stop("Invalid ABIF directory count")
  dirs <- vector("list", n_dir)
  for (i in seq_len(n_dir)) {
    off <- dir_off + (i-1L)*28L
    dirs[[i]] <- parse_dir_entry(raw_slice0(r, off, 28L))
  }
  key <- vapply(dirs, function(e) paste0(e$tag, e$number), character(1))
  names(dirs) <- key
  structure(list(raw=r, version=version, dirs=dirs, path=path), class="abif_file")
}

abif_get_entry <- function(abif, tag, number=1L) {
  abif$dirs[[paste0(tag, as.integer(number))]]
}

abif_get_raw <- function(abif, tag, number=1L, default=NULL) {
  e <- abif_get_entry(abif, tag, number)
  if (is.null(e)) return(default)
  entry_bytes(abif$raw, e)
}

abif_get <- function(abif, tag, number=1L, default=NULL) {
  e <- abif_get_entry(abif, tag, number)
  if (is.null(e)) return(default)
  entry_value(abif$raw, e)
}

read_ab1 <- function(path) {
  a <- read_abif(path)
  seq_raw <- abif_get_raw(a, "PBAS", 2L, NULL)
  if (is.null(seq_raw)) seq_raw <- abif_get_raw(a, "PBAS", 1L, raw())
  # Basecall strings may carry a trailing NUL. Drop NUL bytes before conversion.
  seq_raw <- seq_raw[as.integer(seq_raw) != 0L]
  seq <- if (length(seq_raw)) rawToChar(seq_raw) else ""
  seq <- toupper(gsub("[^ACGTRYKMSWBDHVN]", "N", seq))

  q_raw <- abif_get_raw(a, "PCON", 2L, NULL)
  if (is.null(q_raw)) q_raw <- abif_get_raw(a, "PCON", 1L, raw())
  # PCON is declared as char in ABIF, but its bytes are numeric Phred values,
  # so do not pass through a character encoding.
  q <- as.numeric(as.integer(q_raw))

  ploc <- abif_get(a, "PLOC", 2L, NULL)
  if (is.null(ploc)) ploc <- abif_get(a, "PLOC", 1L, numeric())
  ploc <- as.numeric(ploc)

  order_raw <- abif_get_raw(a, "FWO_", 1L, raw())
  order_raw <- order_raw[as.integer(order_raw) != 0L]
  order <- if (length(order_raw)) rawToChar(order_raw) else ""
  order <- toupper(gsub("[^ACGT]", "", order))
  if (nchar(order) != 4L || !setequal(strsplit(order, "")[[1]], BASES)) stop("Unexpected FWO_1 channel order")

  chans <- setNames(vector("list", 4L), strsplit(order, "")[[1]])
  for (i in 0:3) {
    v <- abif_get(a, "DATA", 9L+i, numeric())
    chans[[substr(order, i+1L, i+1L)]] <- as.numeric(v)
  }

  n <- min(nchar(seq), length(q), length(ploc))
  if (n < 1L) stop("No usable basecalls/qualities/peak locations in AB1")
  seq <- substr(seq, 1L, n)
  q <- q[seq_len(n)]
  ploc <- ploc[seq_len(n)]

  list(
    name=sub("\\.[Aa][Bb]1$", "", basename(path)), path=path,
    seq=seq, qual=q, ploc=ploc, channels=chans, channel_order=order,
    trace_len=max(vapply(chans, length, integer(1))), abif_version=a$version
  )
}
