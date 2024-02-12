For the matlab scripts to work, you need to extract files from
"confounds.tar.gz", eg with "tar xfz *.tar.gz" in linux. This should
produce one *.csv file per participant, with each row in that file 
representing one volume (TR) and the 8 columns  representing x, y, z 
translation, row, pitch, yaw rotations, WM and CSF.
