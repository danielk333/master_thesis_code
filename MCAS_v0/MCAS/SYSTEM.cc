#include <math.h>
#include <iostream>
#include <sys/stat.h>
#include <cstdio>
#include <tr1/memory>

#include "define.hh"
#include "functions.hh"

std::string exec(const char* cmd) {
    std::tr1::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

void remove_folder(const std::string& name) {
  namespace fs = boost::filesystem;
  fs::remove_all(name);
}

std::vector<std::string> list_files(const std::string& folder) {
  std::vector<std::string> ret;
  namespace fs = boost::filesystem;

  boost::progress_timer t( std::clog );

  fs::path full_path( fs::initial_path<fs::path>() );

  full_path = fs::system_complete( fs::path( folder ) );

  unsigned long file_count = 0;
  unsigned long dir_count = 0;
  unsigned long other_count = 0;
  unsigned long err_count = 0;

  if ( !fs::exists( full_path ) ) {
    std::cout << "\nNot found: " << full_path.string() << std::endl;
    return ret;
  }

  if ( fs::is_directory( full_path ) ) {
    std::cout << "\nIn directory: " << full_path.string() << "\n\n";
    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( full_path ); dir_itr != end_iter; ++dir_itr ) {
      try {
        if ( fs::is_directory( dir_itr->status() ) ) {
          ++dir_count;
          std::cout << dir_itr->path().filename() << " [directory]\n";
        }
        else if ( fs::is_regular_file( dir_itr->status() ) ) {
          ++file_count;
          ret.push_back(dir_itr->path().string());
          std::cout << dir_itr->path().filename() << "\n";
        }
        else {
          ++other_count;
          std::cout << dir_itr->path().filename() << " [other]\n";
        }

      }
      catch ( const std::exception & ex ) {
        ++err_count;
        std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
      }
    }
    std::cout << "\n" << file_count << " files\n"
              << dir_count << " directories\n"
              << other_count << " others\n"
              << err_count << " errors\n";
  }
  else { // must be a file
    std::cout << "\nFound: " << full_path.string() << "\n";    
  }
  return ret;
}

bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

bool create_file(std::string name) {

  std::ofstream out_pos(name.c_str(), std::ios::app);
  if(out_pos.fail()) {
    return FALSE;
  }

  out_pos.close();

  return TRUE;
}

bool folder_exists(const std::string& name) {
  struct stat status;   

   if ( access( name.c_str(), 0 ) == 0 ) {
      stat( name.c_str(), &status );

      if ( status.st_mode & S_IFDIR ) {
         return TRUE;
      }
      else {
         return FALSE;
      }
   }
   else {
      return FALSE;
   }
}

std::string get_selfpath() {
    char buff[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      return std::string(buff);
    }
    else {}
    	return "";
    /* handle error condition */
}

