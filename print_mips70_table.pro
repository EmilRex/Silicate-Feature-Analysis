; Written for cchen's reponse to Tushar's paper's review
; By EC on 7/21/14

pro print_mips70_table

; *************************************************** ;

; Load Tushar's fit data from db generated csv
fmt = 's,f,f,f,f,f,f,f,f,f,f,f'
readcol, 'disk_new.csv',F=fmt,disk_name,chisq,rin,rout,rlaw,amin,amax,alaw,mass,fcryst,foliv,ffost

;fmt = ''
;readcol, 'multi_new.csv',F=fmt,multi_name;,params

; Get names associated with MIPS70 values
matches = strip_ext('old_savfiles_mcmc','sav')

; Read in response curve
fmt='f,f,f,f,f,f'
readcol,'MIPS_SED_response.txt',F=fmt,v1,v2,resp_wave,response,v5,v6, /SILENT

; Get Teff, amin and dist_val for the object
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val


; *************************************************** ;
; Loop over each object

FOR i=0, (n_elements(disk_name)-1) DO BEGIN
  FOR j=0, (n_elements(matches)-1) DO BEGIN
  
    IF (disk_name[i] eq matches[j]) THEN BEGIN
  
      ; *************************************************** ;
      ; Load correct data
      
      ; Load MIPS70 data from IDL savefile
      restore, 'old_savfiles_mcmc/'+disk_name[i]+'.sav'
      
      ; Load star system constants
      for k = 0, size(catalog_name,/n_elements)-1 do begin
        if (catalog_name[k] eq disk_name[i]) then begin
          Teff=c_teff[k]
          amin=c_amin[k]
          dist_val=c_dist_val[k]
        endif
      endfor
    
      ; Load global property values
      COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
      COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
      COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit
      restore, 'graintempdata.sav'
      restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene
      
      ; *************************************************** ;
      ; Retrieve grainprops
      
      effectiveTemp = Teff
      
      ; Calculate based on masses based on Isochrones etc
      ; find the right grain model for calculating temperatures
      cmd = 'ls modelgrids/Teff*grains.sav'
      spawn, cmd, grainfiles
      
      ; read in temperatures
      strbeg = strpos(grainfiles, 'Teff')+3
      strend = strpos(grainfiles, 'grains')
      
      tarray = fltarr(n_elements(grainfiles))
      for l=0,n_elements(tarray)-1 do begin
        tarray[l] = float(strmid(grainfiles[l], strbeg[l]+1, strend[l]-strbeg[l]-1))
      endfor
      
      ; reorder arrays
      ii = sort(tarray)
      grainfiles = grainfiles[ii]
      tarray = tarray[ii]
      
      kuruczindex = interpol(findgen(n_elements(tarray)),tarray,effectiveTemp)
      ki = round(kuruczindex) < (n_elements(tarray)-1) > 0
      
      ; next command will retrieve temptable, folivine, Teff
      ; as generated in generateallgrid.pro
      ; folivine = [0.0, 1.0]
      
      restore, grainfiles[ki]
      
      ; *************************************************** ;
      ; Calculate and interpolate modeled spectrum from csv

      ; multi_mips
      ;link[2]=10^(link[2])
      ;link[8]=10^(link[8])
      ;out_model = modeltwograin(model_x, link)
      
      
      ; disk_mips
      ; Convert db data to model data    
      ;link[0]=10^(link[0])
      ;link[6]=10^(link[6])
      diskmass = f(mass)
      ; ...
      
      ; Define fit parameters
      params = [rin,rout,rlaw,amin,amax,alaw,diskmass,foliv,fcryst,ffost]
      
      ; For reference
      ;rin  = params[0]
      ;rout = exp(params[1])*rin
      ;rlaw = params[2]
      ;amin = params[3]
      ;amax = exp(params[4])*amin
      ;alaw = params[5]
      ;diskmass = params[6]
      ;foliv = params[7]
      ;fcrys = params[8]
      ;fforst = params[9]
      
     
      ; Sample model on same wavelengths as response curve
      out_model = modelsinglespectrum(transpose(resp_wave), params)
            
      ; *************************************************** ;
      ; Normalize MIPS SED data by comparing with MIPS70
      ; -> needs some clean-up after it has been tested and works correctly
      
      ; Redefine input names to match names in below
      mips_sed_spec = out_model
      mips_sed_wave = resp_wave
      
      ; Define constants
      c = 3.0e10 ; cm/s
      to_cm = 1.0e-4
      
      ;Define integrand
      int1 = (c*1.0e-23*mips_sed_spec)/(mips_sed_wave*to_cm)^2 ; convert units
      ;int2 = interpol(int1,mips_sed_wave*to_cm,resp_wave*to_cm) ; get flux at right spots
      int2=int1
      int3 = int2*response ; multiply flux and response
      
      ; Find integral of filter
      int_filter = int_tabulated(resp_wave*to_cm,response)
      
      ; Integrate the MIPS SED over the bandpass
      synth_f70 = int_tabulated(resp_wave*to_cm,int3)/int_filter
      
      ; Convert MIPS 70 to proper units
      MIPS_as_flambda = (MIPS70_VAL*c*1.0e-23)/(71.42*to_cm)^2
      
      ; Find the normalization constant
      norm = MIPS_as_flambda/synth_f70
      
      ; Normalize
      mips_sed_spec = norm*mips_sed_spec
      
      ; *************************************************** ;
      ; Print output as latex style table

      ; header: [name?], actual mips 70, actual error, mips 70 according to model

      print,format='(%" %f %f %f ")', $
        MIPS70_VAL, MIPS70_ERR, mips_sed_spec
      ; *************************************************** ;
      
    ENDIF
  ENDFOR
ENDFOR

END