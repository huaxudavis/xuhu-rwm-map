#Script to split haplo file and generate consensus

## Summary ##

This script group spp sequences and generates consensus sequence from loc file.


## Details ##

### Program logic flow ###
  * [split\_flowchart.pdf](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/split_flowchart.pdf)
  * [README\_split\_spps.doc](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/README_split_spps.doc)

### Parameters ###
  * Input loc file
  * window size, window step and max mismatch(optional)

### Sample Command ###
```
$ python split_spps_v4.py haplo_test_50K.loc 0.1 0.05 10

BACP
BACQ
BACR
BACT
BACU
... ...
Please find log file at haplo_test_50K.log. 
```

### Input and Output ###

  * Input file: -- loc file: [haplo\_test\_50K.loc](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/haplo_test_50K.loc)

  * Output file: -- log file: [haplo\_test\_50K.log](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/haplo_test_50K.log)

  * Output group file: -- example file: [BACP\_0.out](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/BACP_0.out)  [BACQ\_0.out](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/BACQ_0.out)   [BACR\_0.out](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/BACR_0.out)  [BACT\_0.out](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/BACT_0.out) [BACU\_0.out](http://xuhu-rwm-map.googlecode.com/svn/wiki/data/BACU_0.out)