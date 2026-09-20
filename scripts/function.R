

check_path <- function(path) {
  sub("/?$", "/", path)
}

add_path <- function(OutPath, filename) {
  paste0(OutPath, filename)
}



