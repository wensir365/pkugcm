# Change Log

## 2017-Apr-9: New Profiling with OpenMP

At pkuclimate.club, a new profiling with 2-thread OpenMP was performed as below. Got any clue? Seems useless...  :(

### Flat profile:

```
Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 22.10     85.43    85.43   843185     0.00     0.00  mktend_
 12.94    135.45    50.02  2659566     0.00     0.00  sp2fc_
 11.35    179.32    43.87   863740     0.00     0.00  dv2uv_
 11.26    222.87    43.54   691165     0.00     0.00  dfft8_
 10.86    264.86    41.99   690938     0.00     0.00  ifft8_
  9.20    300.44    35.59   677265     0.00     0.00  dfft4_
  6.91    327.17    26.73   691045     0.00     0.00  ifft4_
  6.67    352.97    25.80   648842     0.00     0.00  ifft2_
  3.60    366.87    13.90   691167     0.00     0.00  gp2fc_
  2.55    376.73     9.86                             calcgp_
  1.40    382.13     5.40    86403     0.00     0.00  spectral_
  0.58    384.36     2.23    85726     0.00     0.00  sp2fcdmu_
  0.43    386.04     1.68   172806     0.00     0.00  fc2sp_
  0.12    386.50     0.46    86400     0.00     0.00  energy_
  0.02    386.59     0.09    86403     0.00     0.00  gridpoint_
  0.01    386.61     0.02   259209     0.00     0.00  altcs_
  0.00    386.62     0.02                             ifft3_
  0.00    386.63     0.01   642197     0.00     0.00  fc2gp_
  0.00    386.64     0.01   345612     0.00     0.00  mpgallsp_
  0.00    386.65     0.01                             alt2reg_
  0.00    386.66     0.01                             dfft3_
```

### Call graph:

```
granularity: each sample hit covers 2 byte(s) for 0.00% of 386.66 seconds

index % time    self  children    called     name
                               86403             calcgp_ [1]
[1]     92.0    9.86  346.00       0+86403   calcgp_ [1]
               85.43    0.00  843185/843185      mktend_ [3]
               12.16   69.24  604764/691167      gp2fc_ [2]
               50.02    0.00 2659566/2659566     sp2fc_ [5]
                0.01   45.46  555794/642197      fc2gp_ [4]
               43.87    0.00  863740/863740      dv2uv_ [6]
               36.74    0.00  604535/690938      ifft8_ [8]
                2.23    0.00   85726/85726       sp2fcdmu_ [16]
                0.84    0.00   86403/172806      fc2sp_ [17]
                               86403             calcgp_ [1]
-----------------------------------------------
                1.74    9.89   86403/691167      gridpoint_ [14]
               12.16   69.24  604764/691167      calcgp_ [1]
[2]     24.1   13.90   79.13  691167         gp2fc_ [2]
               43.54    0.00  691165/691165      dfft8_ [7]
               35.59    0.00  677265/677265      dfft4_ [9]
-----------------------------------------------
               85.43    0.00  843185/843185      calcgp_ [1]
[3]     22.1   85.43    0.00  843185         mktend_ [3]
-----------------------------------------------
                0.00    7.07   86403/642197      gridpoint_ [14]
                0.01   45.46  555794/642197      calcgp_ [1]
[4]     13.6    0.01   52.53  642197         fc2gp_ [4]
               26.73    0.00  691045/691045      ifft4_ [12]
               25.80    0.00  648842/648842      ifft2_ [13]
-----------------------------------------------
               50.02    0.00 2659566/2659566     calcgp_ [1]
[5]     12.9   50.02    0.00 2659566         sp2fc_ [5]
-----------------------------------------------
               43.87    0.00  863740/863740      calcgp_ [1]
[6]     11.3   43.87    0.00  863740         dv2uv_ [6]
-----------------------------------------------
               43.54    0.00  691165/691165      gp2fc_ [2]
[7]     11.3   43.54    0.00  691165         dfft8_ [7]
-----------------------------------------------
                5.25    0.00   86403/690938      gridpoint_ [14]
               36.74    0.00  604535/690938      calcgp_ [1]
[8]     10.9   41.99    0.00  690938         ifft8_ [8]
-----------------------------------------------
               35.59    0.00  677265/677265      gp2fc_ [2]
[9]      9.2   35.59    0.00  677265         dfft4_ [9]
-----------------------------------------------
                0.00   30.77       1/1           MAIN__ [11]
[10]     8.0    0.00   30.77       1         master_ [10]
                0.09   24.81   86403/86403       gridpoint_ [14]
                5.40    0.01   86403/86403       spectral_ [15]
                0.46    0.00   86400/86400       energy_ [18]
                0.00    0.00    3600/3600        outsp_ [30]
                0.00    0.00    3600/3600        outgp_ [29]
                0.00    0.00    1080/1080        wrzs_ [33]
                0.00    0.00       4/4           makebm_ [45]
-----------------------------------------------
                                                 <spontaneous>
[11]     8.0    0.00   30.77                 MAIN__ [11]
                0.00   30.77       1/1           master_ [10]
                0.00    0.00       1/1           mpstart_ [65]
                0.00    0.00       1/1           resolution_ [74]
                0.00    0.00       1/1           read_resolution_ [71]
                0.00    0.00       1/1           epilog_ [53]
                0.00    0.00       1/1           prolog_ [70]
                0.00    0.00       1/1           allocate_arrays_ [51]
                0.00    0.00       1/1           mpstop_ [66]
-----------------------------------------------
               26.73    0.00  691045/691045      fc2gp_ [4]
[12]     6.9   26.73    0.00  691045         ifft4_ [12]
-----------------------------------------------
               25.80    0.00  648842/648842      fc2gp_ [4]
[13]     6.7   25.80    0.00  648842         ifft2_ [13]
-----------------------------------------------
                0.09   24.81   86403/86403       master_ [10]
[14]     6.4    0.09   24.81   86403         gridpoint_ [14]
                1.74    9.89   86403/691167      gp2fc_ [2]
                0.00    7.07   86403/642197      fc2gp_ [4]
                5.25    0.00   86403/690938      ifft8_ [8]
                0.84    0.00   86403/172806      fc2sp_ [17]
                0.02    0.00  259209/259209      altcs_ [19]
                0.00    0.00  345612/345612      mpsumsc_ [24]
                0.00    0.00  259209/259209      mpgacs_ [25]
                0.00    0.00   86403/86403       mpgagp_ [27]
                0.00    0.00   86403/86403       mpsum_ [28]
-----------------------------------------------
                5.40    0.01   86403/86403       master_ [10]
[15]     1.4    5.40    0.01   86403         spectral_ [15]
                0.01    0.00  345612/345612      mpgallsp_ [21]
-----------------------------------------------
                2.23    0.00   85726/85726       calcgp_ [1]
[16]     0.6    2.23    0.00   85726         sp2fcdmu_ [16]
-----------------------------------------------
                0.84    0.00   86403/172806      gridpoint_ [14]
                0.84    0.00   86403/172806      calcgp_ [1]
[17]     0.4    1.68    0.00  172806         fc2sp_ [17]
-----------------------------------------------
                0.46    0.00   86400/86400       master_ [10]
[18]     0.1    0.46    0.00   86400         energy_ [18]
                0.00    0.00    1440/1440        powerprint_ [31]
-----------------------------------------------
                0.02    0.00  259209/259209      gridpoint_ [14]
[19]     0.0    0.02    0.00  259209         altcs_ [19]
-----------------------------------------------
                                                 <spontaneous>
[20]     0.0    0.02    0.00                 ifft3_ [20]
-----------------------------------------------
```

## 2017-Apr-9: OpenMP works!

With OpenMP directives inserted in the "gridpoint" subroutine of puma.f90, OpenMP paralleling works now! But the performance looks not very scalable, detailed test info on pkuclimate.club are as follows:

Arguments: gfortran with -O3 and -fopenmp options, running for 2 years under T21.

### pkuclimate.club (4 processors)

| gfortran -O3 -fopenmp      | 1 CPU | 2 CPU | 3 CPU | 4 CPU |
|:---------------------------|------:|------:|------:|------:|
| total seconds (2-year/T21) |    99 |    68 |    80 |    80 |
| simulate # yr/day          |  1749 |  2523 |  2180 |  2160 |

Seems 2 CPUs make the best shot! Not satisfied, need to find more room to improve further.

### my desktop PC (4 processors)

| gfortran -O3 -fopenmp      | 1 CPU | 2 CPU | 3 CPU | 4 CPU | 5 CPU | 6 CPU | 7 CPU | 8 CPU |
|:---------------------------|------:|------:|------:|------:|------:|------:|------:|------:|
| total seconds (2-year/T21) |    41 |    28 |    27 |    22 |    34 |    26 |    29 |    27 |
| simulate # yr/day          |  4221 |  6199 |  6460 |  7596 |  5139 |  6552 |  6050 |  6400 |
| total seconds (1-year/T42) |   275 |   185 |   169 |   140 |       |       |       |       |
| simulate # yr/day          |   314 |   467 |   511 |   619 |       |       |       |       |

The speed-up ratio on 4 processors reachs the best of 1.8 (T21) or 1.95 (T42).


## 2017-Apr-8: Clean fftmod.f90, gaussmod.f90, and legsym.f90

Format these 3 files with F95 standard. Try to mark some PURE subroutines in fftmod.f90. I would like to do so to facilitate adding OpenACC directives next.

## 2017-Apr-7: Install PGI Community Edition and test NVIDIA GTX-1060

Install the free community edition of PGI compilers, with support of OpenACC working with NVIDIA GPUs. I received the basic information about my GTX-1060 card via command "pgaccelinfo":

```
wensir@himalaya /work/model/pkugcm/doc $ pgaccelinfo 

CUDA Driver Version:           8000
NVRM version:                  NVIDIA UNIX x86_64 Kernel Module  375.39  Tue Jan 31 20:47:00 PST 2017

Device Number:                 0
Device Name:                   GeForce GTX 1060 3GB
Device Revision Number:        6.1
Global Memory Size:            3153068032
Number of Multiprocessors:     9
Concurrent Copy and Execution: Yes
Total Constant Memory:         65536
Total Shared Memory per Block: 49152
Registers per Block:           65536
Warp Size:                     32
Maximum Threads per Block:     1024
Maximum Block Dimensions:      1024, 1024, 64
Maximum Grid Dimensions:       2147483647 x 65535 x 65535
Maximum Memory Pitch:          2147483647B
Texture Alignment:             512B
Clock Rate:                    1708 MHz
Execution Timeout:             Yes
Integrated Device:             No
Can Map Host Memory:           Yes
Compute Mode:                  default
Concurrent Kernels:            Yes
ECC Enabled:                   No
Memory Clock Rate:             4004 MHz
Memory Bus Width:              192 bits
L2 Cache Size:                 1572864 bytes
Max Threads Per SMP:           2048
Async Engines:                 2
Unified Addressing:            Yes
Managed Memory:                Yes
PGI Compiler Option:           -ta=tesla:cc60
```

## 2017-Apr-4: Add srv-to-nc converter

Add a script "make_netcdf" in ppp/nc to perform converting a service file into netcdf format using Climate Data Operators (CDO). Note that the vertical coordinates, identified as the 2nd integer in the header of a service file, cannot be recognized correctly. Should be "lev", but was converted to "time" instead.

## 2017-Apr-3: Make ppp working as the pre-processor

ppp.f90 works with legsym.f90, fftmod.f90, and gaussmod.f90 in ../src. I re-organized makefile and make it work independantly under current clean framework. Note that use it with following steps:

1. compile the excutable "ppp.x" ("make")
2. prepare namelist "ppp_namelist" ("ln -s namelist/ppp_namelist .")
3. prepare namelist "resolution_namelist" ("ln -s namelist/resolution_namelist .")
4. prepare orography data file "N032_surf_0129.sra" ("ln -s data/N032_surf_0129.sra ." for T21 or "N064_surf_0129.sra" for T42))
5. run the excutable "./ppp.x > log"

Then, you can see 4 addition files as output: 134 for surface pressure; 121 for constant part of temperature profile; 122 for departure of temperature profile; 123 for the relaxiation of rest. temperature. And 1 plain text file ppp-puma.txt, for information transfered to puma in running stage. I currently neglect this file.

Please note that ppp.x could be compatible with T21 and T42, according to the comments in its source code. Not sure if it works fine with other resolutions.

## 2017-Mar-30: Profiling information

For the first time I profiling the code with "-pg" argument in compling phase and using "gprof -b puma.x gmon.out > profile" after performing a 10-yr T21 run. The results are very very helpful, like this:


### Flat profile:

```
Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 15.14     46.87    46.87    86403     0.00     0.00  calcgp_
 14.44     91.59    44.72   864030     0.00     0.00  mktend_
 10.54    124.24    32.65   691224     0.00     0.00  dfft4_
  8.96    151.99    27.75    86403     0.00     0.00  spectral_
  8.92    179.63    27.63   691224     0.00     0.00  dfft8_
  8.30    205.33    25.70   691224     0.00     0.00  ifft8_
  7.45    228.39    23.06   864030     0.00     0.00  dv2uv_
  7.43    251.41    23.02  2678493     0.00     0.00  sp2fc_
  6.29    270.87    19.47   691224     0.00     0.00  ifft4_
  6.22    290.14    19.27   691224     0.00     0.00  ifft2_
  2.55    298.04     7.91   691224     0.00     0.00  gp2fc_
  2.22    304.91     6.87    86403     0.00     0.00  gridpoint_
  0.68    307.03     2.12    86400     0.00     0.00  energy_
  0.45    308.42     1.39   172806     0.00     0.00  fc2sp_
  0.28    309.29     0.87    86403     0.00     0.00  sp2fcdmu_
  0.03    309.38     0.09   259209     0.00     0.00  altcs_
  0.02    309.43     0.05   345612     0.00     0.00  mpgallsp_
  0.02    309.48     0.05   345612     0.00     0.00  mpsumsc_
  0.02    309.53     0.05     3600     0.00     0.00  outsp_
  0.01    309.56     0.04   691224     0.00     0.00  fc2gp_
  0.01    309.58     0.02   147601     0.00     0.00  writesp_
  0.01    309.60     0.02                             dv2uv_alt_
  0.01    309.62     0.02                             filter_spectral_modes_
  0.00    309.64     0.02     1440     0.00     0.00  powerprint_
  0.00    309.65     0.01        1     0.01     0.01  legini_
  0.00    309.66     0.01        1     0.01   309.61  master_
  0.00    309.66     0.01                             ifft3_
  0.00    309.67     0.01                             mpsumval_
```

### Call graph

```
granularity: each sample hit covers 2 byte(s) for 0.00% of 309.67 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]    100.0    0.00  309.62                 MAIN__ [1]
                0.01  309.60       1/1           master_ [2]
                0.00    0.01       1/1           prolog_ [28]
                0.00    0.00       1/1           mpstart_ [67]
                0.00    0.00       1/1           resolution_ [75]
                0.00    0.00       1/1           read_resolution_ [72]
                0.00    0.00       1/1           epilog_ [56]
                0.00    0.00       1/1           allocate_arrays_ [54]
                0.00    0.00       1/1           mpstop_ [68]
-----------------------------------------------
                0.01  309.60       1/1           MAIN__ [1]
[2]    100.0    0.01  309.60       1         master_ [2]
                6.87  272.73   86403/86403       gridpoint_ [3]
               27.75    0.05   86403/86403       spectral_ [9]
                2.12    0.02   86400/86400       energy_ [16]
                0.05    0.02    3600/3600        outsp_ [20]
                0.00    0.00    3600/3600        outgp_ [34]
                0.00    0.00    1080/1080        wrzs_ [36]
                0.00    0.00       4/4           makebm_ [48]
-----------------------------------------------
                6.87  272.73   86403/86403       master_ [2]
[3]     90.3    6.87  272.73   86403         gridpoint_ [3]
                7.91   60.28  691224/691224      gp2fc_ [4]
               46.87    0.00   86403/86403       calcgp_ [5]
               44.72    0.00  864030/864030      mktend_ [6]
                0.04   38.73  691224/691224      fc2gp_ [7]
               25.70    0.00  691224/691224      ifft8_ [11]
               23.06    0.00  864030/864030      dv2uv_ [12]
               23.02    0.00 2678493/2678493     sp2fc_ [13]
                1.39    0.00  172806/172806      fc2sp_ [17]
                0.87    0.00   86403/86403       sp2fcdmu_ [18]
                0.09    0.00  259209/259209      altcs_ [19]
                0.05    0.00  345612/345612      mpsumsc_ [22]
                0.00    0.00  259209/259209      mpgacs_ [31]
                0.00    0.00   86403/86403       mpgagp_ [32]
                0.00    0.00   86403/86403       mpsum_ [33]
-----------------------------------------------
                7.91   60.28  691224/691224      gridpoint_ [3]
[4]     22.0    7.91   60.28  691224         gp2fc_ [4]
               32.65    0.00  691224/691224      dfft4_ [8]
               27.63    0.00  691224/691224      dfft8_ [10]
-----------------------------------------------
               46.87    0.00   86403/86403       gridpoint_ [3]
[5]     15.1   46.87    0.00   86403         calcgp_ [5]
-----------------------------------------------
               44.72    0.00  864030/864030      gridpoint_ [3]
[6]     14.4   44.72    0.00  864030         mktend_ [6]
-----------------------------------------------
                0.04   38.73  691224/691224      gridpoint_ [3]
[7]     12.5    0.04   38.73  691224         fc2gp_ [7]
               19.47    0.00  691224/691224      ifft4_ [14]
               19.27    0.00  691224/691224      ifft2_ [15]
                0.00    0.00       1/1           fftini_ [57]
-----------------------------------------------
               32.65    0.00  691224/691224      gp2fc_ [4]
[8]     10.5   32.65    0.00  691224         dfft4_ [8]
-----------------------------------------------
               27.75    0.05   86403/86403       master_ [2]
[9]      9.0   27.75    0.05   86403         spectral_ [9]
                0.05    0.00  345612/345612      mpgallsp_ [21]
-----------------------------------------------
               27.63    0.00  691224/691224      gp2fc_ [4]
[10]     8.9   27.63    0.00  691224         dfft8_ [10]
-----------------------------------------------
               25.70    0.00  691224/691224      gridpoint_ [3]
[11]     8.3   25.70    0.00  691224         ifft8_ [11]
-----------------------------------------------
               23.06    0.00  864030/864030      gridpoint_ [3]
[12]     7.4   23.06    0.00  864030         dv2uv_ [12]
-----------------------------------------------
               23.02    0.00 2678493/2678493     gridpoint_ [3]
[13]     7.4   23.02    0.00 2678493         sp2fc_ [13]
-----------------------------------------------
               19.47    0.00  691224/691224      fc2gp_ [7]
[14]     6.3   19.47    0.00  691224         ifft4_ [14]
-----------------------------------------------
               19.27    0.00  691224/691224      fc2gp_ [7]
[15]     6.2   19.27    0.00  691224         ifft2_ [15]
-----------------------------------------------
                2.12    0.02   86400/86400       master_ [2]
[16]     0.7    2.12    0.02   86400         energy_ [16]
                0.02    0.00    1440/1440        powerprint_ [26]
-----------------------------------------------
                1.39    0.00  172806/172806      gridpoint_ [3]
[17]     0.4    1.39    0.00  172806         fc2sp_ [17]
-----------------------------------------------
                0.87    0.00   86403/86403       gridpoint_ [3]
[18]     0.3    0.87    0.00   86403         sp2fcdmu_ [18]
-----------------------------------------------
                0.09    0.00  259209/259209      gridpoint_ [3]
[19]     0.0    0.09    0.00  259209         altcs_ [19]
-----------------------------------------------
                0.05    0.02    3600/3600        master_ [2]
[20]     0.0    0.05    0.02    3600         outsp_ [20]
                0.02    0.00  147601/147601      writesp_ [23]
-----------------------------------------------
                0.05    0.00  345612/345612      spectral_ [9]
[21]     0.0    0.05    0.00  345612         mpgallsp_ [21]
-----------------------------------------------
                0.05    0.00  345612/345612      gridpoint_ [3]
[22]     0.0    0.05    0.00  345612         mpsumsc_ [22]
-----------------------------------------------
                0.02    0.00  147601/147601      outsp_ [20]
[23]     0.0    0.02    0.00  147601         writesp_ [23]
-----------------------------------------------
                                                 <spontaneous>
[24]     0.0    0.02    0.00                 dv2uv_alt_ [24]
-----------------------------------------------
                                                 <spontaneous>
[25]     0.0    0.02    0.00                 filter_spectral_modes_ [25]
-----------------------------------------------
                0.02    0.00    1440/1440        energy_ [16]
[26]     0.0    0.02    0.00    1440         powerprint_ [26]
-----------------------------------------------
                0.01    0.00       1/1           prolog_ [28]
[27]     0.0    0.01    0.00       1         legini_ [27]
-----------------------------------------------
                0.00    0.01       1/1           MAIN__ [1]
[28]     0.0    0.00    0.01       1         prolog_ [28]
                0.01    0.00       1/1           legini_ [27]
                0.00    0.00      32/32          mpbcr_ [41]
                0.00    0.00      28/30          mpbci_ [42]
                0.00    0.00      18/18          mpbcrn_ [43]
                0.00    0.00       7/8           mpscsp_ [46]
                0.00    0.00       3/3           mpscrn_ [51]
                0.00    0.00       3/8           mpbcl_ [45]
                0.00    0.00       2/2           mpscdn_ [53]
                0.00    0.00       2/2           mpbcin_ [52]
                0.00    0.00       1/1           initruido_ [63]
                0.00    0.00       1/1           mpscin_ [66]
                0.00    0.00       1/1           initfd_ [60]
                0.00    0.00       1/1           printseed_ [71]
                0.00    0.00       1/1           restart_ini_ [76]
                0.00    0.00       1/1           inilat_ [59]
                0.00    0.00       1/1           inigau_ [58]
                0.00    0.00       1/1           readnl_ [74]
                0.00    0.00       1/1           legpri_ [65]
                0.00    0.00       1/1           initsi_ [64]
                0.00    0.00       1/1           initpm_ [61]
                0.00    0.00       1/1           altlat_ [55]
                0.00    0.00       1/1           initrandom_ [62]
-----------------------------------------------

Index by function name

  [54] allocate_arrays_       [64] initsi_                [30] mpsumval_
  [19] altcs_                 [27] legini_                [69] noise_
  [55] altlat_                [65] legpri_                [35] ntodat_
   [5] calcgp_                [37] lubksb_                [34] outgp_
   [8] dfft4_                 [39] ludcmp_                [20] outsp_
  [10] dfft8_                 [48] makebm_                [26] powerprint_
  [12] dv2uv_                  [2] master_                [70] printprofile_
  [24] dv2uv_alt_             [40] minvers_               [71] printseed_
  [16] energy_                 [6] mktend_                [28] prolog_
  [56] epilog_                [42] mpbci_                 [44] put_restart_array_
   [7] fc2gp_                 [52] mpbcin_                [47] put_restart_integer_
  [17] fc2sp_                 [45] mpbcl_                 [38] ql_
  [57] fftini_                [41] mpbcr_                 [72] read_resolution_
  [25] filter_spectral_modes_ [43] mpbcrn_                [50] read_surf_
   [4] gp2fc_                 [31] mpgacs_                [73] read_vargp_
   [3] gridpoint_             [32] mpgagp_                [74] readnl_
  [15] ifft2_                 [21] mpgallsp_              [75] resolution_
  [29] ifft3_                 [49] mpputsp_               [76] restart_ini_
  [14] ifft4_                 [53] mpscdn_                [77] restart_prepare_
  [11] ifft8_                 [66] mpscin_                [78] set_vertical_grid_
  [58] inigau_                [51] mpscrn_                [79] setzt_
  [59] inilat_                [46] mpscsp_                [13] sp2fc_
  [60] initfd_                [67] mpstart_               [18] sp2fcdmu_
  [61] initpm_                [68] mpstop_                 [9] spectral_
  [62] initrandom_            [33] mpsum_                 [23] writesp_
  [63] initruido_             [22] mpsumsc_               [36] wrzs_
```

## 2017-Mar-29: Make "doc" directory and produce PDF for all sources

Creat "doc" directory and move LICENSE and readme.md into it. Make "make_srcpdf" script to produce a PDF book for all sources in ../src and readme.md.

## 2017-Mar-28: Rename directory puma to src

Rename "puma", the directory of all fortran source code, into "src", to make the directory structure neater.

## 2017-Mar-27: Add MPI function

Add mpimod.f90 and make minor modification on makefile, so that I can deliver MPI runs using 4 CPUs. Seems the excutable with MPI still can be run with 1 single CPU, like this: "./puma_mpi.x 32 10". Normally, the MPI run can be done by "mpiexec -np 4 puma_mpi.x 32 10".

The running speed, under the unif of "simulate model years per day",  was tested on my Intel i5-6400 CPU (roughly 100 GFlops) with 4 processors, with gfortran as the compilor and -O3 as the optimal argument.


| PUMA           |   T21 |  T42 |
|:---------------|------:|-----:|
| 1 CPU (No MPI) |  4373 |  299 |
| 4 CPU (MPI)    | 13918 |  782 |

| PlaSim         |   T21 |  T42 |
|:---------------|------:|-----:|
| 1 CPU (No MPI) |   155 |   23 |
| 4 CPU (MPI)    |   514 |   67 |

Note:

- T21 has 64(lon)x32(lat)=2048 grids; T42 has 128(lon)x64(lat)=8192 grids.
- the 1st number here (4373) is 1755 tested on pkuclimate.club, an 5-year-old Intel i3-530 CPU.

### CPU info for my PC

```
wensir@himalaya /proc $ lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                4
On-line CPU(s) list:   0-3
Thread(s) per core:    1
Core(s) per socket:    4
Socket(s):             1
NUMA node(s):          1
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 94
Model name:            Intel(R) Core(TM) i5-6400 CPU @ 2.70GHz
Stepping:              3
CPU MHz:               3099.937
CPU max MHz:           3300.0000
CPU min MHz:           800.0000
BogoMIPS:              5424.13
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              6144K
NUMA node0 CPU(s):     0-3
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush 
dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs 
bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx est 
tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave 
avx f16c rdrand lahf_lm abm 3dnowprefetch intel_pt tpr_shadow vnmi flexpriority ept vpid fsgsbase 
tsc_adjust bmi1 avx2 smep bmi2 erms invpcid mpx rdseed adx smap clflushopt xsaveopt xsavec xgetbv1 
dtherm ida arat pln pts hwp hwp_notify hwp_act_window hwp_epp
```

### CPU info for climateserver.3322.org (pkuclimate.club)

```
wensir@pkuclimate:~$ lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                4
On-line CPU(s) list:   0-3
Thread(s) per core:    2
Core(s) per socket:    2
Socket(s):             1
NUMA node(s):          1
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 37
Model name:            Intel(R) Core(TM) i3 CPU         530  @ 2.93GHz
Stepping:              2
CPU MHz:               2926.231
BogoMIPS:              5852.46
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              4096K
NUMA node0 CPU(s):     0-3
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush 
dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx rdtscp lm constant_tsc arch_perfmon pebs bts rep_good 
nopl xtopology nonstop_tsc aperfmperf pni dtes64 monitor ds_cpl vmx est tm2 ssse3 cx16 xtpr pdcm 
sse4_1 sse4_2 popcnt lahf_lm tpr_shadow vnmi flexpriority ept vpid dtherm arat
```

## 2017-Mar-25: Add Post-Processing (pp)

Add directory "pp" for the post-processing excutable "burn7.x" and 2 associated namelist, one for 3 surface variables, another for 10 multi-level variables. To compile burn7.x you should apt install libnetcdf-cxx-legacy-dev first.

## 2017-Mar-25: Remove GUI-related code
Remove guimod_stub.f90, and comment out all the GUI-related code in puma.f90 with a description "XW(Mar/25/2017) to remove GUI:". Please note I did not delete any lines, just comment out. The number of source code were shinked from 7 to 6.

## 2017-Mar-23: Make PUMA running

Successfully make PUMA running on 1CPU, with fake MPI and GUI interface (stub). 
I changed "nresources" line in subroutine "epilog" of puma.f90, which is related to retrieving process time, memory, and disk info used in pumax.c. They are just printed into puma_diag file at the closing procedure of a run. No big deal. The code are marked with "XW".

## 2017-Mar-17: Make PLASIM17 running with MOST

### PlaSim 4core 1year T21 run

```
wensir@himalaya /work/model/puma17 $ ./most.x 
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  resmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  plasimmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  radmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  oceanmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  icemod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  seamod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  cpl_stub.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod_stub.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  rainmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  landmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  simba.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  p_earth.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  plasim.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  calmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  outmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  miscmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fluxmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  surfmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  tracermod.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  tpcore.f90
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  trc_routines.f90
mpicc -c -O3 pumax_stub.c
mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  lsgmod.f90
mpif90 -o plasim.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan mpimod.o 
fftmod.o guimod_stub.o rainmod.o simba.o p_earth.o resmod.o plasim.o plasimmod.o calmod.o gaussmod.o 
legmod.o outmod.o miscmod.o fluxmod.o radmod.o surfmod.o landmod.o seamod.o icemod.o oceanmod.o 
restartmod.o tracermod.o tpcore.o trc_routines.o pumax_stub.o lsgmod.o cpl_stub.o 

=== Success: Launched process most_plasim_run ===
```

### PUMA 4core 1year T21 run

	wensir@himalaya ~/model/puma17 $ ./most.x 
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  puma.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod.f90
	mpicc -c -O3 pumax.c
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legini.f90
	as -o legfast32.o legfast32.s
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
	mpif90 -o puma.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan 
               mpimod.o fftmod.o guimod.o pumax.o legini.o legfast32.o puma.o restartmod.o gaussmod.o 
               -L/usr/lib/X11 -lX11

	=== Success: Launched process most_puma_run ===

### PlaSim 1core 1year T21 run

```
wensir@himalaya /work/model/puma17 $ ./most.x 
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  resmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  plasimmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod_stub.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  radmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  oceanmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  icemod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  seamod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  cpl_stub.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod_stub.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  rainmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  landmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  simba.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  p_earth.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  plasim.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  calmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  outmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  miscmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fluxmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  surfmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  tracermod.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  tpcore.f90
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  trc_routines.f90
gcc   -c -O3 -I /usr/lib/X11/include pumax_stub.c
gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  lsgmod.f90
gfortran -o plasim.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan mpimod_stub.o 
fftmod.o guimod_stub.o rainmod.o simba.o p_earth.o resmod.o plasim.o plasimmod.o calmod.o gaussmod.o 
legmod.o outmod.o miscmod.o fluxmod.o radmod.o surfmod.o landmod.o seamod.o icemod.o oceanmod.o 
restartmod.o tracermod.o tpcore.o trc_routines.o pumax_stub.o lsgmod.o cpl_stub.o 

=== Success: Launched process most_plasim_run ===
```

### PUMA 1core 1year T21 run

	wensir@himalaya ~/model/puma17 $ ./most.x 
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  puma.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod_stub.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod.f90
	gcc   -c -O3 -I /usr/lib/X11/include pumax.c
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legini.f90
	as -o legfast32.o legfast32.s
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
	gfortran -o puma.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan 
                 mpimod_stub.o fftmod.o guimod.o pumax.o legini.o legfast32.o puma.o restartmod.o gaussmod.o 
                 -L/usr/lib/X11 -lX11

	=== Success: Launched process most_puma_run ===


### Configure and Compile Source Code

	wensir@himalaya ~/model/puma17 $ ./configure.sh 
	Found FORTRAN-90 compiler     at: /usr/bin/gfortran
	gfortran version  5.4
	Found Xlib (X11)              at: /usr/lib/X11
	Found C compiler              at: /usr/bin/gcc
	Found C++ compiler            at: /usr/bin/c++
	Found MPI FORTRAN-90 compiler at: /usr/bin/mpif90
	Found MPI C compiler          at: /usr/bin/mpicc
	Found MPI C++ compiler        at: /usr/bin/mpicxx
	Found GNU assembler           at: /usr/bin/as
	Fast Legendre Transformation    : ACTIVE (Linux)
	gfortran -o f90check.x f90check.f90
	gcc -o cc_check.x cc_check.c

	System info for <himalaya>
	Architecture: Linux himalaya 4.4.0-53-generic #74-Ubuntu SMP Fri Dec 2 15:59:10 UTC 2016 x86_64 GNU/Linux
	Endian format             : little endian
	FORTRAN control word size : 4 bytes
	FORTRAN integer size      : 4 bytes
	FORTRAN real    size      : 4 bytes
	C       int     size      : 4 bytes
	C       float   size      : 4 bytes
	C       long    size      : 8 bytes
	C  long long    size      : 8 bytes
	FORTRAN Compiler: gfortran
	C       Compiler: gcc
	gcc -o most.x most.c -I/usr/lib/X11/include -lm -L/usr/lib/X11 -lX11

	configuration complete - run <most.x>


### Install required library following README_UBUNTU

- apt install libx11-dev
- apt install openmpi-bin openmpi-common openmpi-doc libopenmpi-dev
