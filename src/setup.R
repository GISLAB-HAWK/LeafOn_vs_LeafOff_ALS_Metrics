#-------------------------------------------------------------
# Name:         setup.R
# Description:  Script sets up a working environment,
#               defines file paths for data import and output,
#               and loads required packages.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------



# 01 - setup working environment
#--------------------------------

# create directory called 'data' with sub directories
# 'raw_data', 'processed_data', and 'metadata'
if (!file.exists(paste('data')) |
    (!file.exists(paste('data/raw_data')) |
     (!file.exists(paste('data/raw_data/forest_inventory')) |
      (!file.exists(paste('data/raw_data/pc_leafoff_2024')) |
       (!file.exists(paste('data/raw_data/pc_leafon_2023')) ||
        (!file.exists(paste('data/processed_data')) |
         (!file.exists(paste('data/processed_data/forest_inventory')) |
          (!file.exists(paste('data/processed_data/chm_leafoff')) |
           (!file.exists(paste('data/processed_data/chm_leafon')) |
            (!file.exists(paste('data/processed_data/gap_polygons_leafoff')) |
             (!file.exists(paste('data/processed_data/gap_polygons_leafon')) |
              (!file.exists(paste('data/processed_data/metrics')) |
               (!file.exists(paste('data/processed_data/models')) |
                (!file.exists(paste('data/processed_data/predictions')) |
                 (!file.exists(paste('data/metadata'))
                  )))))))))))))))
  {
  
  dir.create('data')
  dir.create('data/raw_data')
  dir.create('data/raw_data/forest_inventory')
  dir.create('data/raw_data/pc_leafoff_2024')
  dir.create('data/raw_data/pc_leafon_2023')
  dir.create('data/processed_data')
  dir.create('data/processed_data/forest_inventory')
  dir.create('data/processed_data/chm_leafoff')
  dir.create('data/processed_data/chm_leafon')
  dir.create('data/processed_data/gap_polygons_leafoff')
  dir.create('data/processed_data/gap_polygons_leafon')
  dir.create('data/processed_data/metrics')
  dir.create('data/processed_data/models')
  dir.create('data/processed_data/predictions')
  dir.create('data/metadata')
  
} else {
  
  invisible()
  
}

# create directory called 'src'
if (!file.exists(paste('src'))) {
  
  dir.create('src')
  
} else {
  
  invisible()
  
}

# create directory called 'docs'
if (!file.exists(paste('docs'))) {
  
  dir.create('docs')
  
} else {
  
  invisible()
  
}

# create directory called 'scripts'
if (!file.exists(paste('scripts'))) {
  
  dir.create('scripts')
  
} else {
  
  invisible()
  
}

# create directory called 'output'
if (!file.exists(paste('output'))) {
  
  dir.create('output')
  
} else {
  
  invisible()
  
}

# list the files and directories
list.files(recursive = TRUE, include.dirs = TRUE)



# 02 - file path definitions
#---------------------------

# define raw data directory
raw_data_dir <- 'data/raw_data/'

# define processed data directory
processed_data_dir <- 'data/processed_data/'

# define meta data directory
meta_data_dir <- 'data/metadata/'

# define output directory
output_dir <- 'output/'



# 03 - package loading
#----------------------

# load (and install) required packages
load_packages <- function(packages, github_remotes = NULL, github_repos = NULL, gitlab_remotes = NULL) {
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Package '", pkg, "' not found, attempting to install from CRAN...", sep = ""))
      install.packages(pkg, dependencies = TRUE)
      
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop(paste("Package '", pkg, "' not found and could not be installed from CRAN.", sep = ""))
      }
    }
  }
  
  if (!is.null(github_remotes)) {
    for (pkg_name in names(github_remotes)) {
      if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Package '", pkg_name, "' not found, attempting to install from GitHub using remotes (", github_remotes[[pkg_name]], ")...", sep = ""))
        remotes::install_github(github_remotes[[pkg_name]])
        
        if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
          stop(paste("Package '", pkg_name, "' not found and could not be installed from GitHub using remotes.", sep = ""))
        }
      }
    }
  }
  
  if (!is.null(github_repos)) {
    for (pkg_name in names(github_repos)) {
      if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Package '", pkg_name, "' not found, attempting to install from GitHub repository (", github_repos[[pkg_name]], ")...", sep = ""))
        install.packages(pkg_name, repos = github_repos[[pkg_name]])
        
        if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
          stop(paste("Package '", pkg_name, "' not found and could not be installed from GitHub repository.", sep = ""))
        }
      }
    }
  }
  
  if (!is.null(gitlab_remotes)) {
    for (pkg_name in names(gitlab_remotes)) {
      if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Package '", pkg_name, "' not found, attempting to install from GitLab using remotes (", gitlab_remotes[[pkg_name]]$repo, ")...", sep = ""))
        
        repo <- gitlab_remotes[[pkg_name]]$repo
        build_vignettes <- if (!is.null(gitlab_remotes[[pkg_name]]$build_vignettes)) {
          gitlab_remotes[[pkg_name]]$build_vignettes
        } else {
          FALSE
        }
        
        remotes::install_gitlab(repo, build_vignettes = build_vignettes)
        
        if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
          stop(paste("Package '", pkg_name, "' not found and could not be installed from GitLab.", sep = ""))
        }
      }
    }
  }
}

load_packages(
  c('terra', 'lidR' , 'sf', 'stats','dplyr', 'ggplot2',
    'mgcv', 'scam', 'cowplot', 'ggrepel', 'caret', 'CAST',
    'parallel', 'doParallel', 'ranger'),
  github_remotes = c(TreeGrOSSinR = 'rnuske/TreeGrOSSinR'),
  github_repos = c(lasR = 'https://r-lidar.r-universe.dev'),
  gitlab_remotes = list(
    rBDAT = list(repo = 'vochr/rBDAT', build_vignettes = TRUE)
  )
)