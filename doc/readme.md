# Change Log

## 2017-Mar-30: Profiling information

For the first time I profiling the code with "-pg" argument in compling phase and using "gprof -b puma.x gmon.out > profile" after performing a 10-yr T21 run. The results are very very helpful, like this:


### Flat profile:

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
  0.00    309.67     0.00   259209     0.00     0.00  mpgacs_
  0.00    309.67     0.00    86403     0.00     0.00  mpgagp_
  0.00    309.67     0.00    86403     0.00     0.00  mpsum_
  0.00    309.67     0.00     3600     0.00     0.00  outgp_
  0.00    309.67     0.00     1080     0.00     0.00  ntodat_
  0.00    309.67     0.00     1080     0.00     0.00  wrzs_
  0.00    309.67     0.00      840     0.00     0.00  lubksb_
  0.00    309.67     0.00      310     0.00     0.00  ql_
  0.00    309.67     0.00       84     0.00     0.00  ludcmp_
  0.00    309.67     0.00       84     0.00     0.00  minvers_
  0.00    309.67     0.00       32     0.00     0.00  mpbcr_
  0.00    309.67     0.00       30     0.00     0.00  mpbci_
  0.00    309.67     0.00       18     0.00     0.00  mpbcrn_
  0.00    309.67     0.00       13     0.00     0.00  put_restart_array_
  0.00    309.67     0.00        8     0.00     0.00  mpbcl_
  0.00    309.67     0.00        8     0.00     0.00  mpscsp_
  0.00    309.67     0.00        5     0.00     0.00  put_restart_integer_
  0.00    309.67     0.00        4     0.00     0.00  makebm_
  0.00    309.67     0.00        4     0.00     0.00  mpputsp_
  0.00    309.67     0.00        4     0.00     0.00  read_surf_
  0.00    309.67     0.00        3     0.00     0.00  mpscrn_
  0.00    309.67     0.00        2     0.00     0.00  mpbcin_
  0.00    309.67     0.00        2     0.00     0.00  mpscdn_
  0.00    309.67     0.00        1     0.00     0.00  allocate_arrays_
  0.00    309.67     0.00        1     0.00     0.00  altlat_
  0.00    309.67     0.00        1     0.00     0.00  epilog_
  0.00    309.67     0.00        1     0.00     0.00  fftini_
  0.00    309.67     0.00        1     0.00     0.00  inigau_
  0.00    309.67     0.00        1     0.00     0.00  inilat_
  0.00    309.67     0.00        1     0.00     0.00  initfd_
  0.00    309.67     0.00        1     0.00     0.00  initpm_
  0.00    309.67     0.00        1     0.00     0.00  initrandom_
  0.00    309.67     0.00        1     0.00     0.00  initruido_
  0.00    309.67     0.00        1     0.00     0.00  initsi_
  0.00    309.67     0.00        1     0.00     0.00  legpri_
  0.00    309.67     0.00        1     0.00     0.00  mpscin_
  0.00    309.67     0.00        1     0.00     0.00  mpstart_
  0.00    309.67     0.00        1     0.00     0.00  mpstop_
  0.00    309.67     0.00        1     0.00     0.00  noise_
  0.00    309.67     0.00        1     0.00     0.00  printprofile_
  0.00    309.67     0.00        1     0.00     0.00  printseed_
  0.00    309.67     0.00        1     0.00     0.01  prolog_
  0.00    309.67     0.00        1     0.00     0.00  read_resolution_
  0.00    309.67     0.00        1     0.00     0.00  read_vargp_
  0.00    309.67     0.00        1     0.00     0.00  readnl_
  0.00    309.67     0.00        1     0.00     0.00  resolution_
  0.00    309.67     0.00        1     0.00     0.00  restart_ini_
  0.00    309.67     0.00        1     0.00     0.00  restart_prepare_
  0.00    309.67     0.00        1     0.00     0.00  set_vertical_grid_
  0.00    309.67     0.00        1     0.00     0.00  setzt_

### Call graph

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
                                                 <spontaneous>
[29]     0.0    0.01    0.00                 ifft3_ [29]
-----------------------------------------------
                                                 <spontaneous>
[30]     0.0    0.01    0.00                 mpsumval_ [30]
-----------------------------------------------
                0.00    0.00  259209/259209      gridpoint_ [3]
[31]     0.0    0.00    0.00  259209         mpgacs_ [31]
-----------------------------------------------
                0.00    0.00   86403/86403       gridpoint_ [3]
[32]     0.0    0.00    0.00   86403         mpgagp_ [32]
-----------------------------------------------
                0.00    0.00   86403/86403       gridpoint_ [3]
[33]     0.0    0.00    0.00   86403         mpsum_ [33]
-----------------------------------------------
                0.00    0.00    3600/3600        master_ [2]
[34]     0.0    0.00    0.00    3600         outgp_ [34]
-----------------------------------------------
                0.00    0.00    1080/1080        wrzs_ [36]
[35]     0.0    0.00    0.00    1080         ntodat_ [35]
-----------------------------------------------
                0.00    0.00    1080/1080        master_ [2]
[36]     0.0    0.00    0.00    1080         wrzs_ [36]
                0.00    0.00    1080/1080        ntodat_ [35]
-----------------------------------------------
                0.00    0.00     840/840         minvers_ [40]
[37]     0.0    0.00    0.00     840         lubksb_ [37]
-----------------------------------------------
                0.00    0.00     310/310         inigau_ [58]
[38]     0.0    0.00    0.00     310         ql_ [38]
-----------------------------------------------
                0.00    0.00      84/84          minvers_ [40]
[39]     0.0    0.00    0.00      84         ludcmp_ [39]
-----------------------------------------------
                0.00    0.00      84/84          makebm_ [48]
[40]     0.0    0.00    0.00      84         minvers_ [40]
                0.00    0.00     840/840         lubksb_ [37]
                0.00    0.00      84/84          ludcmp_ [39]
-----------------------------------------------
                0.00    0.00      32/32          prolog_ [28]
[41]     0.0    0.00    0.00      32         mpbcr_ [41]
-----------------------------------------------
                0.00    0.00       2/30          read_resolution_ [72]
                0.00    0.00      28/30          prolog_ [28]
[42]     0.0    0.00    0.00      30         mpbci_ [42]
-----------------------------------------------
                0.00    0.00      18/18          prolog_ [28]
[43]     0.0    0.00    0.00      18         mpbcrn_ [43]
-----------------------------------------------
                0.00    0.00      13/13          epilog_ [56]
[44]     0.0    0.00    0.00      13         put_restart_array_ [44]
-----------------------------------------------
                0.00    0.00       1/8           read_vargp_ [73]
                0.00    0.00       3/8           prolog_ [28]
                0.00    0.00       4/8           read_surf_ [50]
[45]     0.0    0.00    0.00       8         mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       1/8           initfd_ [60]
                0.00    0.00       7/8           prolog_ [28]
[46]     0.0    0.00    0.00       8         mpscsp_ [46]
-----------------------------------------------
                0.00    0.00       5/5           epilog_ [56]
[47]     0.0    0.00    0.00       5         put_restart_integer_ [47]
-----------------------------------------------
                0.00    0.00       4/4           master_ [2]
[48]     0.0    0.00    0.00       4         makebm_ [48]
                0.00    0.00      84/84          minvers_ [40]
-----------------------------------------------
                0.00    0.00       4/4           epilog_ [56]
[49]     0.0    0.00    0.00       4         mpputsp_ [49]
-----------------------------------------------
                0.00    0.00       4/4           initfd_ [60]
[50]     0.0    0.00    0.00       4         read_surf_ [50]
                0.00    0.00       4/8           mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       3/3           prolog_ [28]
[51]     0.0    0.00    0.00       3         mpscrn_ [51]
-----------------------------------------------
                0.00    0.00       2/2           prolog_ [28]
[52]     0.0    0.00    0.00       2         mpbcin_ [52]
-----------------------------------------------
                0.00    0.00       2/2           prolog_ [28]
[53]     0.0    0.00    0.00       2         mpscdn_ [53]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[54]     0.0    0.00    0.00       1         allocate_arrays_ [54]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[55]     0.0    0.00    0.00       1         altlat_ [55]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[56]     0.0    0.00    0.00       1         epilog_ [56]
                0.00    0.00      13/13          put_restart_array_ [44]
                0.00    0.00       5/5           put_restart_integer_ [47]
                0.00    0.00       4/4           mpputsp_ [49]
                0.00    0.00       1/1           restart_prepare_ [77]
-----------------------------------------------
                0.00    0.00       1/1           fc2gp_ [7]
[57]     0.0    0.00    0.00       1         fftini_ [57]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[58]     0.0    0.00    0.00       1         inigau_ [58]
                0.00    0.00     310/310         ql_ [38]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[59]     0.0    0.00    0.00       1         inilat_ [59]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[60]     0.0    0.00    0.00       1         initfd_ [60]
                0.00    0.00       4/4           read_surf_ [50]
                0.00    0.00       1/1           read_vargp_ [73]
                0.00    0.00       1/1           setzt_ [79]
                0.00    0.00       1/1           printprofile_ [70]
                0.00    0.00       1/8           mpscsp_ [46]
                0.00    0.00       1/1           noise_ [69]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[61]     0.0    0.00    0.00       1         initpm_ [61]
                0.00    0.00       1/1           set_vertical_grid_ [78]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[62]     0.0    0.00    0.00       1         initrandom_ [62]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[63]     0.0    0.00    0.00       1         initruido_ [63]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[64]     0.0    0.00    0.00       1         initsi_ [64]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[65]     0.0    0.00    0.00       1         legpri_ [65]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[66]     0.0    0.00    0.00       1         mpscin_ [66]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[67]     0.0    0.00    0.00       1         mpstart_ [67]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[68]     0.0    0.00    0.00       1         mpstop_ [68]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[69]     0.0    0.00    0.00       1         noise_ [69]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[70]     0.0    0.00    0.00       1         printprofile_ [70]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[71]     0.0    0.00    0.00       1         printseed_ [71]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[72]     0.0    0.00    0.00       1         read_resolution_ [72]
                0.00    0.00       2/30          mpbci_ [42]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[73]     0.0    0.00    0.00       1         read_vargp_ [73]
                0.00    0.00       1/8           mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[74]     0.0    0.00    0.00       1         readnl_ [74]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[75]     0.0    0.00    0.00       1         resolution_ [75]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[76]     0.0    0.00    0.00       1         restart_ini_ [76]
-----------------------------------------------
                0.00    0.00       1/1           epilog_ [56]
[77]     0.0    0.00    0.00       1         restart_prepare_ [77]
-----------------------------------------------
                0.00    0.00       1/1           initpm_ [61]
[78]     0.0    0.00    0.00       1         set_vertical_grid_ [78]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[79]     0.0    0.00    0.00       1         setzt_ [79]
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

## 2017-Mar-29: Make "doc" directory and produce PDF for all sources

Creat "doc" directory and move LICENSE and readme.md into it. Make "make_srcpdf" script to produce a PDF book for all sources in ../src and readme.md.

## 2017-Mar-28: Rename directory puma to src

Rename "puma", the directory of all fortran source code, into "src", to make the directory structure neater.

## 2017-Mar-27: Add MPI function

Add mpimod.f90 and make minor modification on makefile, so that I can deliver MPI runs using 4 CPUs. Seems the excutable with MPI still can be run with 1 single CPU, like this: "./puma_mpi.x 32 10". Normally, the MPI run can be done by "mpiexec -np 4 puma_mpi.x 32 10".

The running speed, under the unif of "simulate model years per day",  was tested on my Intel i5 CPU with 4 processors, with gfortran as the compilor and -O3 as the optimal argument.

| `gfortran -O3` |   T21 |  T42 |
|:---------------|------:|-----:|
| 1 CPU (No MPI) |  4373 |  299 |
| 4 CPU (MPI)    | 13918 |  782 |

Note that the 1st number here (4373) is 1755 tested on pkuclimate.club, an 5-year-old Intel i3 Xeon CPU.

## 2017-Mar-25: Add Post-Processing (pp)

Add directory "pp" for the post-processing excutable "burn7.x" and 2 associated namelist, one for 3 surface variables, another for 10 multi-level variables. To compile burn7.x you should apt install libnetcdf-cxx-legacy-dev first.

## 2017-Mar-25: Remove GUI-related code
Remove guimod_stub.f90, and comment out all the GUI-related code in puma.f90 with a description "XW(Mar/25/2017) to remove GUI:". Please note I did not delete any lines, just comment out. The number of source code were shinked from 7 to 6.

## 2017-Mar-23: Make PUMA running

Successfully make PUMA running on 1CPU, with fake MPI and GUI interface (stub). 
I changed "nresources" line in subroutine "epilog" of puma.f90, which is related to retrieving process time, memory, and disk info used in pumax.c. They are just printed into puma_diag file at the closing procedure of a run. No big deal. The code are marked with "XW".

## 2017-Mar-17: Make PLASIM17 running with MOST

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

Flat profile:

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
  0.00    309.67     0.00   259209     0.00     0.00  mpgacs_
  0.00    309.67     0.00    86403     0.00     0.00  mpgagp_
  0.00    309.67     0.00    86403     0.00     0.00  mpsum_
  0.00    309.67     0.00     3600     0.00     0.00  outgp_
  0.00    309.67     0.00     1080     0.00     0.00  ntodat_
  0.00    309.67     0.00     1080     0.00     0.00  wrzs_
  0.00    309.67     0.00      840     0.00     0.00  lubksb_
  0.00    309.67     0.00      310     0.00     0.00  ql_
  0.00    309.67     0.00       84     0.00     0.00  ludcmp_
  0.00    309.67     0.00       84     0.00     0.00  minvers_
  0.00    309.67     0.00       32     0.00     0.00  mpbcr_
  0.00    309.67     0.00       30     0.00     0.00  mpbci_
  0.00    309.67     0.00       18     0.00     0.00  mpbcrn_
  0.00    309.67     0.00       13     0.00     0.00  put_restart_array_
  0.00    309.67     0.00        8     0.00     0.00  mpbcl_
  0.00    309.67     0.00        8     0.00     0.00  mpscsp_
  0.00    309.67     0.00        5     0.00     0.00  put_restart_integer_
  0.00    309.67     0.00        4     0.00     0.00  makebm_
  0.00    309.67     0.00        4     0.00     0.00  mpputsp_
  0.00    309.67     0.00        4     0.00     0.00  read_surf_
  0.00    309.67     0.00        3     0.00     0.00  mpscrn_
  0.00    309.67     0.00        2     0.00     0.00  mpbcin_
  0.00    309.67     0.00        2     0.00     0.00  mpscdn_
  0.00    309.67     0.00        1     0.00     0.00  allocate_arrays_
  0.00    309.67     0.00        1     0.00     0.00  altlat_
  0.00    309.67     0.00        1     0.00     0.00  epilog_
  0.00    309.67     0.00        1     0.00     0.00  fftini_
  0.00    309.67     0.00        1     0.00     0.00  inigau_
  0.00    309.67     0.00        1     0.00     0.00  inilat_
  0.00    309.67     0.00        1     0.00     0.00  initfd_
  0.00    309.67     0.00        1     0.00     0.00  initpm_
  0.00    309.67     0.00        1     0.00     0.00  initrandom_
  0.00    309.67     0.00        1     0.00     0.00  initruido_
  0.00    309.67     0.00        1     0.00     0.00  initsi_
  0.00    309.67     0.00        1     0.00     0.00  legpri_
  0.00    309.67     0.00        1     0.00     0.00  mpscin_
  0.00    309.67     0.00        1     0.00     0.00  mpstart_
  0.00    309.67     0.00        1     0.00     0.00  mpstop_
  0.00    309.67     0.00        1     0.00     0.00  noise_
  0.00    309.67     0.00        1     0.00     0.00  printprofile_
  0.00    309.67     0.00        1     0.00     0.00  printseed_
  0.00    309.67     0.00        1     0.00     0.01  prolog_
  0.00    309.67     0.00        1     0.00     0.00  read_resolution_
  0.00    309.67     0.00        1     0.00     0.00  read_vargp_
  0.00    309.67     0.00        1     0.00     0.00  readnl_
  0.00    309.67     0.00        1     0.00     0.00  resolution_
  0.00    309.67     0.00        1     0.00     0.00  restart_ini_
  0.00    309.67     0.00        1     0.00     0.00  restart_prepare_
  0.00    309.67     0.00        1     0.00     0.00  set_vertical_grid_
  0.00    309.67     0.00        1     0.00     0.00  setzt_

			Call graph


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
                                                 <spontaneous>
[29]     0.0    0.01    0.00                 ifft3_ [29]
-----------------------------------------------
                                                 <spontaneous>
[30]     0.0    0.01    0.00                 mpsumval_ [30]
-----------------------------------------------
                0.00    0.00  259209/259209      gridpoint_ [3]
[31]     0.0    0.00    0.00  259209         mpgacs_ [31]
-----------------------------------------------
                0.00    0.00   86403/86403       gridpoint_ [3]
[32]     0.0    0.00    0.00   86403         mpgagp_ [32]
-----------------------------------------------
                0.00    0.00   86403/86403       gridpoint_ [3]
[33]     0.0    0.00    0.00   86403         mpsum_ [33]
-----------------------------------------------
                0.00    0.00    3600/3600        master_ [2]
[34]     0.0    0.00    0.00    3600         outgp_ [34]
-----------------------------------------------
                0.00    0.00    1080/1080        wrzs_ [36]
[35]     0.0    0.00    0.00    1080         ntodat_ [35]
-----------------------------------------------
                0.00    0.00    1080/1080        master_ [2]
[36]     0.0    0.00    0.00    1080         wrzs_ [36]
                0.00    0.00    1080/1080        ntodat_ [35]
-----------------------------------------------
                0.00    0.00     840/840         minvers_ [40]
[37]     0.0    0.00    0.00     840         lubksb_ [37]
-----------------------------------------------
                0.00    0.00     310/310         inigau_ [58]
[38]     0.0    0.00    0.00     310         ql_ [38]
-----------------------------------------------
                0.00    0.00      84/84          minvers_ [40]
[39]     0.0    0.00    0.00      84         ludcmp_ [39]
-----------------------------------------------
                0.00    0.00      84/84          makebm_ [48]
[40]     0.0    0.00    0.00      84         minvers_ [40]
                0.00    0.00     840/840         lubksb_ [37]
                0.00    0.00      84/84          ludcmp_ [39]
-----------------------------------------------
                0.00    0.00      32/32          prolog_ [28]
[41]     0.0    0.00    0.00      32         mpbcr_ [41]
-----------------------------------------------
                0.00    0.00       2/30          read_resolution_ [72]
                0.00    0.00      28/30          prolog_ [28]
[42]     0.0    0.00    0.00      30         mpbci_ [42]
-----------------------------------------------
                0.00    0.00      18/18          prolog_ [28]
[43]     0.0    0.00    0.00      18         mpbcrn_ [43]
-----------------------------------------------
                0.00    0.00      13/13          epilog_ [56]
[44]     0.0    0.00    0.00      13         put_restart_array_ [44]
-----------------------------------------------
                0.00    0.00       1/8           read_vargp_ [73]
                0.00    0.00       3/8           prolog_ [28]
                0.00    0.00       4/8           read_surf_ [50]
[45]     0.0    0.00    0.00       8         mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       1/8           initfd_ [60]
                0.00    0.00       7/8           prolog_ [28]
[46]     0.0    0.00    0.00       8         mpscsp_ [46]
-----------------------------------------------
                0.00    0.00       5/5           epilog_ [56]
[47]     0.0    0.00    0.00       5         put_restart_integer_ [47]
-----------------------------------------------
                0.00    0.00       4/4           master_ [2]
[48]     0.0    0.00    0.00       4         makebm_ [48]
                0.00    0.00      84/84          minvers_ [40]
-----------------------------------------------
                0.00    0.00       4/4           epilog_ [56]
[49]     0.0    0.00    0.00       4         mpputsp_ [49]
-----------------------------------------------
                0.00    0.00       4/4           initfd_ [60]
[50]     0.0    0.00    0.00       4         read_surf_ [50]
                0.00    0.00       4/8           mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       3/3           prolog_ [28]
[51]     0.0    0.00    0.00       3         mpscrn_ [51]
-----------------------------------------------
                0.00    0.00       2/2           prolog_ [28]
[52]     0.0    0.00    0.00       2         mpbcin_ [52]
-----------------------------------------------
                0.00    0.00       2/2           prolog_ [28]
[53]     0.0    0.00    0.00       2         mpscdn_ [53]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[54]     0.0    0.00    0.00       1         allocate_arrays_ [54]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[55]     0.0    0.00    0.00       1         altlat_ [55]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[56]     0.0    0.00    0.00       1         epilog_ [56]
                0.00    0.00      13/13          put_restart_array_ [44]
                0.00    0.00       5/5           put_restart_integer_ [47]
                0.00    0.00       4/4           mpputsp_ [49]
                0.00    0.00       1/1           restart_prepare_ [77]
-----------------------------------------------
                0.00    0.00       1/1           fc2gp_ [7]
[57]     0.0    0.00    0.00       1         fftini_ [57]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[58]     0.0    0.00    0.00       1         inigau_ [58]
                0.00    0.00     310/310         ql_ [38]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[59]     0.0    0.00    0.00       1         inilat_ [59]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[60]     0.0    0.00    0.00       1         initfd_ [60]
                0.00    0.00       4/4           read_surf_ [50]
                0.00    0.00       1/1           read_vargp_ [73]
                0.00    0.00       1/1           setzt_ [79]
                0.00    0.00       1/1           printprofile_ [70]
                0.00    0.00       1/8           mpscsp_ [46]
                0.00    0.00       1/1           noise_ [69]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[61]     0.0    0.00    0.00       1         initpm_ [61]
                0.00    0.00       1/1           set_vertical_grid_ [78]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[62]     0.0    0.00    0.00       1         initrandom_ [62]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[63]     0.0    0.00    0.00       1         initruido_ [63]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[64]     0.0    0.00    0.00       1         initsi_ [64]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[65]     0.0    0.00    0.00       1         legpri_ [65]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[66]     0.0    0.00    0.00       1         mpscin_ [66]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[67]     0.0    0.00    0.00       1         mpstart_ [67]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[68]     0.0    0.00    0.00       1         mpstop_ [68]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[69]     0.0    0.00    0.00       1         noise_ [69]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[70]     0.0    0.00    0.00       1         printprofile_ [70]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[71]     0.0    0.00    0.00       1         printseed_ [71]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[72]     0.0    0.00    0.00       1         read_resolution_ [72]
                0.00    0.00       2/30          mpbci_ [42]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[73]     0.0    0.00    0.00       1         read_vargp_ [73]
                0.00    0.00       1/8           mpbcl_ [45]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[74]     0.0    0.00    0.00       1         readnl_ [74]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[75]     0.0    0.00    0.00       1         resolution_ [75]
-----------------------------------------------
                0.00    0.00       1/1           prolog_ [28]
[76]     0.0    0.00    0.00       1         restart_ini_ [76]
-----------------------------------------------
                0.00    0.00       1/1           epilog_ [56]
[77]     0.0    0.00    0.00       1         restart_prepare_ [77]
-----------------------------------------------
                0.00    0.00       1/1           initpm_ [61]
[78]     0.0    0.00    0.00       1         set_vertical_grid_ [78]
-----------------------------------------------
                0.00    0.00       1/1           initfd_ [60]
[79]     0.0    0.00    0.00       1         setzt_ [79]
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
