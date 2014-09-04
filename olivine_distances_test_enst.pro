PRO olivine_distances_test_enst
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

; Load names and distances from db produced csv
fmt = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'
readcol, 'multi_new_test_enst.csv',F=fmt,db_name,chisq,temp1,temp2,Loc1,Loc2,amin1,amin2,mass1,mass2,fcryst1,fcryst2,oliv1,oliv2,ffost1,ffost2

agr = transpose([[amin1],[amin2]])
temp = transpose([[temp1],[temp2]])
foliv = transpose([[oliv1],[oliv2]])
fcrys = transpose([[fcryst1],[fcryst2]])
ffors = transpose([[ffost1],[ffost2]])
loc = transpose([[Loc1],[Loc2]])

; Load star system parameters: Teff and dist
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val, /silent

; Create save file and header
close,1
openw,1,'olivine_distances_test_enst.csv'
printf,1, 'Name,Temp1,Tushar_Dist1,Oliv_Dist1,Pyro_Dist1,Fors_Dist1,Enst_Dist1,Temp2,Tushar_Dist2,Oliv_Dist2,Pyro_Dist2,Fors_Dist2,Enst_Dist2'

FOR i=0, n_elements(db_name)-1 DO BEGIN
  
  ; ******************************** ;
  ; Get Data
  
  dist = dblarr(2,4)
  
  FOR part=0,1 DO BEGIN
  
    FOR type=0,3 DO BEGIN
      
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
      qlookup_old, [agr[part,i]], final_phot_wave, foliv[part,i], fcrys[part,i], ffors[part,i], qabs_phot, /separate
      
      ; ******************************** ;
      ; Calculate
    
      int_const = (c*(dist_val_AU^2))/(((final_phot_wave*mm_to_cm)^2))
      
      numerator = INT_TABULATED(final_phot_wave*mm_to_cm,qabs_phot[*,*,type]*int_const*final_phot_fnu*1.0e-23)
    
      blambda = ( (2.0*h*(c^2))/((final_phot_wave*mm_to_cm)^5) )/( exp( (h*c)/((final_phot_wave*mm_to_cm)*k*temp[part,i]) ) -1.0 )
      denominator = 4.0*INT_TABULATED((final_phot_wave*mm_to_cm),qabs_phot[*,*,type]*blambda)
        
      dist[part,type] = sqrt( numerator/denominator )
    
      ENDFOR ; grain type
    
    ENDFOR ; cold,hot disk
    
    printf,1,string(db_name[i])+','+string(temp[0,i])+','+string(loc[0,i])+','+string(dist[0,0])+','+string(dist[0,1])+','+string(dist[0,2])+$
      ','+string(dist[0,3])+','+string(temp[1,i])+','+string(loc[1,i])+','+string(dist[1,0])+','+$
      string(dist[1,1])+','+string(dist[1,2])+','+string(dist[1,3])
      
    ;printf,1,string(db_name[i])+','+string(temp[0,i])+','+string(dist)+','+string(loc[part,i])+','+string(dist-loc)
    
    print, string(db_name[i])+' Done!'

    ;stop
  ENDFOR ; object

close,1

END
