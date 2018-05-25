NULL

.onLoad <- function(libname, pkgname) {

  osname <- get_os()
  ini_file <- get_ini_path()
  morphr_settings <- ini::read.ini(ini_file)
  
  if (osname == "Darwin") 
    morphr_settings$path$dp_shape_match_path <- system.file('bin', 'DPShapeMatch', package = 'morphr')

  if (osname == "unix")
    morphr_settings$path$dp_shape_match_path <- system.file('bin', 'DPShapeMatch_ubuntu_x86_64', package = 'morphr')
  
  ini::write.ini(morphr_settings, ini_file)  

}
