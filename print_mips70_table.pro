; Written for cchen's reponse to Tushar's paper's review
; By EC on 7/21/14

pro print_mips70_table, disk=disk, multi=multi

; *************************************************** ;
; Load Tushar's fit data from db generated csv

if keyword_set(disk) then begin
  out_type = 'disk'
  fmt = 'a,f,f,f,f,f,f,f,f,f,f,f'
  readcol, 'disk_new.csv',F=fmt,name,chisq,rin,rout,rlaw,amin,amax,alaw,diskmass,fcryst,foliv,ffost
endif

if keyword_set(multi) then begin
  out_type = 'multi'
  fmt = 'a,a,f,f,f,f,f,f,f,f,f,f,f,f,f'
  readcol, 'multi.csv',F=fmt,name,type,chisq,temp1,amin1,mass1,oliv1,fcryst1,ffost1,temp2,amin2,mass2,oliv2,fcryst2,ffost2
endif


; *************************************************** ;
; Load Stellar parameters

; Load R-star values
fmt = 'a,f'
readcol, 'r_star.csv',F=fmt,r_star_name,c_r_star

; Get names associated with MIPS70 values
matches = strip_ext('old_savfiles_mcmc','sav')

; Read in response curve
fmt='f,f,f,f,f,f'
readcol,'MIPS_SED_response.txt',F=fmt,v1,v2,resp_wave,response,v5,v6, /SILENT

; Get Teff, amin and dist_val for the object
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val

; Open output file
close,1
openw,1,'print_mips70_table_'+out_type+'.csv'

; Print header
;out_format='(%"%s & %s & %s & %s \\\\")'
out_format='(%"%s, %s, %s, %s, %s, %s")'
printf,1,format=out_format, $
  'Name', 'Chisq', 'Observed MIPS70', 'Observed MIPS70 Error', 'Modeled MIPS70', 'Interpolated MIPS70'

; Define constants
c = 3.0e10 ; cm/s
to_cm = 1.0e-4
m_moon = 7.34767309e22 ; in kg, from google
r_sun = 0.00464913034 ;AU
au_in_cm = 1.49597871e13
pc_in_cm = 3.08567758e18

; *************************************************** ;
; Loop over each object

FOR i=0, (n_elements(name)-1) DO BEGIN
  FOR j=0, (n_elements(matches)-1) DO BEGIN
  
    IF (name[i] eq matches[j]) THEN BEGIN
  
      ; *************************************************** ;
      ; Load correct data
      
      ; Load MIPS70 data from IDL savefile
      restore, 'old_savfiles_mcmc/'+name[i]+'.sav'
      
      ; Load star system constants
      for k = 0, size(catalog_name,/n_elements)-1 do begin
        if (catalog_name[k] eq name[i]) then begin
          Teff=c_teff[k]
          dist_val=c_dist_val[k]
        endif
      endfor

      ; Load correct r_star
      for k=0, n_elements(r_star_name)-1 do begin
        if (r_star_name[k] eq name[i]) then begin
          r_star = c_r_star[k]
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
      
      
      if keyword_set(multi) then begin
        ; *************************************************** ;
        ; Calculate multi spectrum from csv
  
        link = [temp1[i],amin1[i],mass1[i],fcryst1[i],oliv1[i],ffost1[i],$
               temp2[i],amin2[i],mass2[i],fcryst2[i],oliv2[i],ffost2[i]]
        
        ; Sample model on same wavelengths as response curve
        ; Outputs flux in Jansky's        
        out_model = modeltwograin_old(resp_wave, link)
      
      endif
      
      if keyword_set(disk) then begin

        ; *************************************************** ;
        ; Calculate disk spectrum from csv
        
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
        
        ; Convert db values to fit parameters
        r_star_in_au = r_star*r_sun ; r_star in AU
        fit_rin = alog10(rin[i]/(r_sun*r_star)) ;for HD3003 ;rin[i]
        fit_rout = alog(rout[i]/rin[i])
        fit_amax = alog(amax[i]/amin[i])
        
        fit_diskmass = alog10((diskmass[i])*m_moon) ; convert m/m_moon to log(mass_in_kg)
        ; need to apply correction Tushar mentioned
        
        ; Combine fit parameters
        params = [10^fit_rin,fit_rout,rlaw[i],amin[i],fit_amax,alaw[i],10^fit_diskmass,foliv[i],fcryst[i],ffost[i]]
        
        ; Sample model on same wavelengths as response curve
        ; Outputs flux in Jansky's
        out_model = modelsinglespectrum(transpose(resp_wave), params)

      endif

      ; *************************************************** ;
      ; Normalize MIPS SED data by comparing with MIPS70

      ;Define integrand
      int=out_model*response
      
      ; Find integral of filter
      int_filter = int_tabulated(resp_wave*to_cm,response)
      
      ; Integrate the MIPS SED over the bandpass
      synth_f70 = int_tabulated(resp_wave*to_cm,int)/int_filter
      ;synth_f70 = interpol(out_model,resp_wave*to_cm,71.42*to_cm)
      
      interpol_f70 = interpol(out_model,resp_wave*to_cm,71.42*to_cm)
      
      ; *************************************************** ;
      ; Print output as latex style table

      ; header: name, actual mips 70, actual error, mips 70 according to model
      ;out_format='(%"%s & %s & %s & %s \\\\")'
      out_format='(%"%s, %f, %f, %f, %f, %f")'
      printf,1,format=out_format, $
      name[i], chisq[i], MIPS70_VAL, MIPS70_ERROR, synth_f70, interpol_f70
      
      print, string((i+1))+': '+string(name[i])+' done'
      ; *************************************************** ;
      
    ENDIF
  ENDFOR
ENDFOR

close,1

END