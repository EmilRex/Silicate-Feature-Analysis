; +
; NAME
;  pso_min_example
;
; PURPOSE
;  Search the parameter space with the PSO routine rmd_pso 
;
; INPUTS
;  AMIN:
;  TEFF: 
;  NAME:
;  OUT_PAR:
;  SEQUENCE:
;
; KEYWORDS
;   NONE
;
; OUTPUTS
;  OUTPUT: array containing likelihood (i.e. fitness) values 
;
; AUTHORS
;  Tushar Mittal 
;  Christine Chen
;  Emil Christensen - chris2er@dukes.jmu.edu
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY
;  Adapted from rmd_pso.pro by TM (June 2013) as pso_min_example_m.pro
;  Organized and commented by EC (6/27/2014)
;  Generalized and renamed by EC (7/14/14)
; -
; *************************************************** ;
; 
; 
; ****************************************************************************************************** ;
function f1_eval,p,_EXTRA = extra
; ****************************************************************************************************** ;
; 
; This function has a couple of minima and maxima.  This
; particular function is the one that the algorithm will
; call and evaluate agents positions, one at a time.

; Create global variables relating to silicate features
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, Qwaterice, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, waterice_emit, effectiveTempArray, stellar_emit
COMMON file_path, in_dir, out_dir, fit_name, object_name

; *************************************************** ;
; MODEL: ONE GRAIN

IF (fit_name eq 'single') THEN BEGIN
  
  link=[p[0],p[1],p[2],p[3],p[4],p[5],p[6]]
  spectra = modelsinglespectrum(transpose(extra.wave),link, /single )
  
ENDIF

; *************************************************** ;
; MODEL: TWO GRAIN

IF (fit_name eq 'multi') THEN BEGIN 

  link=[p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13]]
  spectra = modelsinglespectrum(transpose(extra.wave),link, /multi )                         

ENDIF

; *************************************************** ;
; MODEL: CONTINUOUS DISK

IF (fit_name eq 'disk') THEN BEGIN

  link=[p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10]]
  spectra = modelsinglespectrum(transpose(extra.wave),link, /disk ) 
  
ENDIF

; *************************************************** ;
; Compute chisq

chisq = TOTAL ( ((extra.spec-spectra)^2.0)/((.05*extra.error)^2.0+(extra.error)^2.0))
z = (chisq)/(extra.dof - 1.)   ; Log of likelihood func - where likelihood funct is exp(-chisq/2 )


; Print to txt file if desired
;openw, 1, 'output_v1/'+extra.name+'_fit_multi_pso.txt',/append
;printf,1,'param,',strtrim(string(p[0]),2),',',strtrim(string(p[1]),2),',',strtrim(string(p[2]),2),',',strtrim(string(p[3]),2),',',strt;rim(string(p[4]),2),',',strtrim(string(p[5]),2),',',strtrim(string(p[6]),2),',',strtrim(string(p[7]),2),',',strtrim(string(p[8]),2),',',strtrim(string(p[9]),2),',',strtrim(string(p[10]),2),',',strtrim(string(p[11]),2),',',strtrim(string(z),2)
;close,1

; Return likelihood and exit
return,z
end


; ****************************************************************************************************** ;
pro pso_test_iterproc,  func, p, iter, interrupt, functargs = functargs, oref = opso, _Extra = extra
; ****************************************************************************************************** ;

COMMON file_path, in_dir, out_dir, fit_name, object_name

compile_opt hidden,idl2
opso->get_property,fresult=fresult

; Write current PSO result to txt file
openw, 1, out_dir+'/'+functargs.name+'_fit_'+fit_name+'_pso_best_'+strtrim(fix(functargs.sequence),2)+'.txt',/append
writeu,1,'Iteration: '+strtrim(string(iter),2)
writeu,1,'fresult : ',fresult

; Write current parameter values to same txt file
for i = 0,n_elements(p)-1 do begin
  strout = 'p['+strtrim(string(i),2)+']='+strtrim(string(p[i]),2)
  writeu,1,strout
endfor

close,1 ; close file
end


; ****************************************************************************************************** ;
pro pso_min_example,amin=amin,Teff=teff_val,name=name1,output=output1,out_par=functargs,sequence=sequence
; ****************************************************************************************************** ;

COMMON file_path, in_dir, out_dir, fit_name, object_name

; This branch "forced_dist" forces ring radius to 90.5 AU -> intended for HD181327
; param = alog10(dist/(r_sun*r_star))
; rstar = 0.60794045 rstar ; rsun = 0.00464913034 AU
; 90.5 AU -> 4.50542

IF (fit_name eq 'single') THEN BEGIN
  ;prange = [[30.0, 800.0],[amin, 30.0],[16.5, 23.5],[0, 1.0],[0, 1.0],[0, 1.0]]
  prange = [[4.50542, 4.50542],[amin, 30.0],[16.5, 23.5],[0, 1.0],[0, 1.0],[0, 1.0],[0, 1.0]]
ENDIF

IF (fit_name eq 'multi') THEN BEGIN
  ;prange = [[30.0, 300.0],[amin, 30.0],[16.5, 23.5],[0., 1.0],[0., 1.0],[0., 1.0],$
  ;         [100.0, 1000.0],[amin, 30.0],[16.5, 23.5],[0., 1.0],[0., 1.0],[0., 1.0]]
  
  prange = [[4.50542, 4.50542],[amin, 30.0],[16.5, 23.5],[0., 1.0],[0., 1.0],[0., 1.0],[0, 1.0],$
           [1.0, 6.0],[amin, 30.0],[16.5, 23.5],[0., 1.0],[0., 1.0],[0., 1.0],[0, 1.0]]
ENDIF

IF (fit_name eq 'disk') THEN BEGIN
  prange = [[1.0, 5.0],[0.0, 10.0],[-5., 5.],[amin, 30.0],[0.0, 7.0],[1.5, 5.5],[16.5, 23.5],[0.0, 1.0],[0.0, 1.0],[0.0, 1.0],[0, 1.0]]
ENDIF

func = 'f1_eval' ; Function to be minimized
n = 40 ; Number of agents in the swarm

; *************************************************** ;
; Store data 

restore,in_dir+'/'+name1+'.sav'

err_chk=where( final_specerr le .01*final_spec)
final_specerr(err_chk)=.01*final_spec[err_chk]
;data_base=[transpose(final_wave),transpose(final_spec),transpose(final_specerr)]
dof=n_elements(final_wave)-n_elements(prange(0,*)); + 1.0  ; Added 1 to account for the extra degree due to mips70

functargs =  {  wave:final_wave,    $
                spec:final_spec,    $
                error:final_specerr,$
                dof:dof,            $
                name:name1,         $
                sequence: sequence }

; Store desired outcome (i.e. actual spectrum)
 data_base=[transpose(final_wave),transpose(final_spec),transpose(final_specerr)]
 fxhmake,header1,data_base,/date
 header_lines = strarr(2)
 header_lines[0] = '/  PSO Fit -  Dust Model - Jang Condell et al. 2013, Mittal et al. 2013'
 header_lines[1] = '/  Spitzer IRS spectrum, Chen et al. 2013'
 sxaddhist, header_lines,header1
 file=out_dir+'/'+name1+'_chn_pso_'+fit_name+'_part_'+strtrim(fix(sequence),1)+'.fits'
 FITS_WRITE,file,data_base,header1
 undefine,data_base;

 ; *************************************************** ;

 ; Create global variables relating to silicate features
 COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, Qwaterice, crystallineabs
 COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
 COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, waterice_emit, effectiveTempArray, stellar_emit
 
 ; Fill the above global variables
 restore, 'graintempdata.sav'
 restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene
 restore, 'qwaterice.sav'

 rho_s = 3.3
 AU_in_cm = 1.496e13

; *************************************************** ;
; Restore silicate features data                                              

; Find the right grain model for calculating temperatures
Teff=teff_val
effectiveTemp = Teff

cmd = 'ls modelgrids/Teff*grains.sav'
spawn, cmd, grainfiles

; read in temperatures                                                                                                     
strbeg = strpos(grainfiles, 'Teff')+3
strend = strpos(grainfiles, 'grains')
tarray = fltarr(n_elements(grainfiles))
for i=0,n_elements(tarray)-1 do begin
   tarray[i] = float(strmid(grainfiles[i], strbeg[i]+1, strend[i]-strbeg[i]-1))
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


; Call the PSO minimization routine
p = rmd_pso(       ftol = .8,                             $
                     function_name = func,                  $
                     FUNCTARGS = functargs,                 $
                     function_value = fval,                 $
                     ncalls = ncalls,                       $
                     weights = [2.051,2.051],               $
                     itmax = 200,                           $
                     quiet = 0B,                            $
                     iterproc = 'pso_test_iterproc',        $
                     vel_fraction = .5,                     $
                     num_particles = n,                     $
                     vel_decrement = .729,                  $
                     prange = prange                        )

print, ''
print, systime()
print,'Result: ',p
print,'Value: ',(fval)
print,'# function calls: ',ncalls

; Rename output and terminate program
output1=[p,fval]
return
end