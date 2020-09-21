REM LUMPREM output files are deleted.
 
del lr_abc.out
del lr_def.out
del lr_ghi.out
 
REM LUMPREM models are run.
 
lumprem lr_abc.in lr_abc.out
lumprem lr_def.in lr_def.out
lumprem lr_ghi.in lr_ghi.out
