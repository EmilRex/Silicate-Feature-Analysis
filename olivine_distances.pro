PRO olivine_distances
; Find the distance given the temperature for all objects in table 2 (olivine only)

; Define constants
pc_in_AU = 2.0626e5
AU_in_cm = 1.496e13
mm_to_cm = 1.0e-4
c = 2.998d10 ; cm/s
h = 6.626d-27 ; erg s
k = 1.381d-16 ; erg/K

; Set graintype
; 0 - Olivine
; 1 - Pyroxine
; 2 - Forsterite
; 3 - Enstatite
type = 0

; Load names and distances from db produced csv
fmt = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'
readcol, 'multi_new.csv',F=fmt,db_name,chisq,temp1,temp2,Loc1,Loc2,amin1,amin2,mass1,mass2,fcryst1,fcryst2,oliv1,oliv2,ffost1,ffost2

; Load star system parameters: Teff and dist
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val, /silent

; Create save file and header
close,1
openw,1,'olivine_distances.csv'
printf,1, 'Name,Temp,Dist,Loc,Diff'

FOR i=0, n_elements(db_name)-1 DO BEGIN
  
  ; ******************************** ;
  ; Get Data
  
  agr = amin1[i]
  temp=temp1[i]
  foliv = oliv1[i]
  fcrys = fcryst1[i]
  ffors = ffost1[i]
  loc = Loc1[i]
  
  ; Load stellar photosphere model data
  restore, 'old_savfiles_mcmc/'+db_name[i]+'.sav'
  
  for j = 0, n_elements(catalog_name)-1 do begin
    if (catalog_name[j] eq db_name[i]) then begin
      Teff=c_teff[j]
      dist_val=c_dist_val[j]
    endif
  endfor
  
  ; Convert distance in parsecs to au
  dist_val_AU = dist_val*pc_in_AU
  
  ; Get q_abs over phot lambdas
  qlookup, [agr], final_phot_wave, foliv, fcrys, ffors, qabs_phot, /separate
  
  ; ******************************** ;
  ; Calculate

  int_const = (c*(dist_val_AU^2))/(((final_phot_wave*mm_to_cm)^2))
  
  numerator = INT_TABULATED(final_phot_wave*mm_to_cm,qabs_phot[*,*,type]*int_const*final_phot_fnu*1.0e-23)

  blambda = ( (2.0*h*(c^2))/((final_phot_wave*mm_to_cm)^5) )/( exp( (h*c)/((final_phot_wave*mm_to_cm)*k*temp) ) -1.0 )
  denominator = 4.0*INT_TABULATED((final_phot_wave*mm_to_cm),qabs_phot[*,*,type]*blambda)
    
  dist = sqrt( numerator/denominator )
  
  printf,1,string(db_name[i])+','+string(temp)+','+string(dist)+','+string(loc)+','+string(dist-loc)
  ;print,string(db_name[i])+' Temp = '+string(temp)+' --- Dist = '+string(dist)
  
ENDFOR

close,1

END