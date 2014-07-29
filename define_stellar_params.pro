pro define_stellar_params

COMMON file_path, in_dir, out_dir, fit_name, object_name
COMMON star_params, t_star, dist_to_star, star_lambda, star_fnu

; Load text file with names, Teff, amin and dist
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val,/silent

; Match loaded data names with name of object being analyzed 
amin = .1
for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq object_name) then begin
       Teff=c_teff[i]
       amin=c_amin[i]
       dist_val=c_dist_val[i]
    endif
endfor

; Load stellar photosphere model data
restore,in_dir+'/'+object_name+'.sav'

t_star = Teff
dist_to_star = dist_val
star_lambda = final_phot_wave
star_fnu = final_phot_fnu

end