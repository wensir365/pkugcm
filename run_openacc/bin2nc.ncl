;==============================================================================================;
; PKU-GCM Post-processing NCL Script                                                           ;
; Developed by Xinyu Wen (xwen@pku.edu.cn), Apr/29/2017                                        ;
; Given key parameters, convert all (years) output files in binary to netCDF for this variable ;
;==============================================================================================;

begin

; ========== 你只需要修改以下关键参数 ==========
;      NLEV, NLON, NTIME, VARID, VARDS, VARUN

   ;--- key parameters (DEMO)
   NLEV  = 1                     ; number of layers
   NLAT  = 32                    ; number of latitudes
   NTIME = 360                   ; default full-length (in days) of a binary file
   VARID = "ps"                  ; variable's name
   VARDS = "Surface Pressure"    ; variable's description
   VARUN = "hPa"                 ; variable's unit

   ;--- key parameters (DEMO)
   NLEV  = 10
   NLAT  = 32
   NTIME = 360
   VARID = "t"
   VARDS = "Temperature"
   VARUN = "K"

   ;--- key parameters (DEMO)
   NLEV  = stringtoint(getenv("NLEV"))
   NLAT  = stringtoint(getenv("NLAT"))
   NTIME = stringtoint(getenv("NTIME"))
   VARID = getenv("VARID")
   VARDS = getenv("VARDS")
   VARUN = getenv("VARUN")

; ========== 你只需要修改以上关键参数 ==========



   ;--- print welcome
   print("                                ")
   print("================================")
   print(" PKU-GCM Post-Processing Script ")
   print(" https://pkugcm.wordpress.com/  ")
   print("--------------------------------")
   print(" VAR Name : "+VARID              )
   print(" VAR Desc : "+VARDS              )
   print(" VAR Unit : "+VARUN              )
   print("--------------------------------")
   print("                                ")

   ;--- constants
   Psurf = 1011.00         ; global mean surface pressure

   ;--- times
   time           = ispan(1,NTIME,1)
   time@long_name = "time"
   time@units     = "days since model starting point 0000-01-00"

   ;--- levels
   if (NLEV.gt.1) then
      sigmh          = fspan(0.0,1.0,NLEV+1)
      sigma          = 0.5*(sigmh(1:NLEV)+sigmh(0:NLEV-1))
      if (VARID.eq."z") then
         lev         = sigmh(0:NLEV-1:-1)
         lev(NLEV-1) = sigma(0)     ; for top layer, only account for half
      else
         lev         = sigma(::-1)
      end if
      lev@long_name  = "sigma coordinate"
      lev@units      = "unitless"
      lev@sigma      = sigma
      lev@sigma_interface = sigmh
   else
      lev            = Psurf
      lev@long_name  = "pressure"
      lev@units      = "hPa"
   end if

   ;--- latitudes
   gau   = tofloat(gaus(NLAT/2))

   lat   = gau(:,0)  ; -90 to 90
   lat@long_name  = "latitude"
   lat@units      = "degrees_north"

   latw  = gau(:,1)
   latw@long_name = "weight for latitude"
   latw@units     = "unitless"

   ;--- longitudes
   NLON     = NLAT*2
   dx       = 360.0/NLON
   dxhalf   = dx/2.0

   lon   = fspan(dxhalf,360.0-dxhalf,NLON)
   lon@long_name  = "longitude"
   lon@units      = "degrees_east"

   ;--- input & output
   filetimestamp  = 0
   do while fileexists("output_"+VARID+"_"+sprinti("%0.4i",filetimestamp)+".bin")
      ; input
      fn    = "output_"+VARID+"_"+sprinti("%0.4i",filetimestamp)+".bin"
      xtmp  = cbinread(fn,-1,"float")
      nt    = dimsizes(xtmp)/(NLEV*NLAT*NLON)
      x     = onedtond(xtmp,(/nt,NLEV,NLAT,NLON/))
      x     = (/ x(:,::-1,::-1,:) /)
      print(fn+": found "+nt+" x ("+NLEV+"x"+NLAT+"x"+NLON+") values")

      ; get var ready
      time  = time + filetimestamp*NTIME

      x!0   = "time"
      x!1   = "lev"
      x!2   = "lat"
      x!3   = "lon"
      x&time= time(0:nt-1)
      x&lev = lev
      x&lat = lat
      x&lon = lon

      x@long_name = VARDS
      x@units     = VARUN

      ; output
      fon   = "output_"+VARID+"_"+sprinti("%0.4i",filetimestamp)+".nc"
      system("rm -f "+fon)
      fo    = addfile(fon,"c")
      filedimdef(fo,"time",-1,True)
      fo->$VARID$ = x
      fo@model    = "Peking University General Circulation Model (PKU-GCM)"
      fo@website  = "find more info via https://pkugcm.wordpress.com/"
      fo@script   = "post-processing script in NCL, written by xwen@pku.edu.cn"

      print(fon+" : generating netcdf file successfully")

      ; loop
      filetimestamp = filetimestamp + 1
      delete(x)
      delete(xtmp)
      print(" ")
   end do

end
