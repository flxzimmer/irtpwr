

# from Remote -------------------------------------------------------------

# Copy files
file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/evaluation_run2.R", "Y:/evaluation_run2.R", overwrite = T)

file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/evaluation_setup2.R", "Y:/evaluation_setup2.R", overwrite = T)


file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/archive/evaluation_setup2_temp2.R", "Y:/evaluation_setup2_temp2.R", overwrite = T)

file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/evaluation_run_temp3.R", "Y:/evaluation_run_temp3.R", overwrite = T)


file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata", "Y:/simsetup.Rdata", overwrite = T)

file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage_0.1.0.tar.gz", "Y:/mmlpwrpackage_0.1.0.tar.gz", overwrite = T)

# from server -------------------------------------------------------------

# install Source package, run
install.packages("mmlpwrpackage_0.1.0.tar.gz",repos=NULL)

nohup nice -n 11 Rscript evaluation_run2.R > run.out &

nohup nice -n 11 Rscript evaluation_setup2.R > run.out &

nohup nice -n 9 Rscript evaluation_setup2_temp2.R > run.out &

nohup nice -n 11 Rscript evaluation_run_temp.R > run.out &

nohup nice -n 11 Rscript evaluation_run_temp3.R > run.out &

pkill -u fzimmer
















# other -------------------------------------------------------------------

# file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage_0.1.0.zip", "Y:/mmlpwrpackage_0.1.0.zip", overwrite = T) # binary

# file.copy("C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/evaluation_run_temp.R", "Y:/evaluation_run_temp.R", overwrite = T)

# nohup nice Rscript evaluation_run_temp.R > run.out &


# devtools::install_local("mmlpwrpackage_0.1.0.zip") # binary
