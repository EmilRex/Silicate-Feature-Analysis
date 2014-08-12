; Separated from main.pro by EC on 7/21/14
; ****************************************************************************************************** ;
function display_historic,name
; ****************************************************************************************************** ;
; Print Tushar's old simulation results to CL

COMMON file_path, in_dir, out_dir, fit_name, object_name

; Set-up default. i.e. missing data
head = "No prior data available"
old_result = " "
chisq_old = " "

; Define constants
m_moon = 7.34767309e25 ; in g, from google
r_sun = 0.00464913034 ;AU

print, ' -------------------------------------------------------------------------------------------- '
print, systime()
print, name,' --- ',fit_name,' --- OLD'
print, ' '

CASE fit_name OF
  'single': BEGIN
    ; No data hence nothing to do
  END
  
  'multi_mips': BEGIN
    fmt = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'
    readcol, 'multi_new.csv',F=fmt,db_name,chisq,temp1,temp2,Loc,Loc2,amin1,amin2,mass1,mass2,fcryst1,fcryst2,oliv1,oliv2,ffost1,ffost2,/silent
    
    FOR i=0,(n_elements(db_name)-1) DO BEGIN
      IF (name eq db_name[i]) THEN BEGIN
        head = ['Temp1 ','a_grain1 ','dustmass1 ','folive1 ','fcrys1 ','ffors1 ','Temp2 ','a_grain2 ','dustmass2 ','folive2 ','fcrys2 ','ffors2 ']
        old_result = [temp1[i],amin1[i],alog10(mass1[i]*m_moon),oliv1[i],fcryst1[i],ffost1[i],temp2[i],amin2[i],alog10(mass2[i]*m_moon),oliv2[i],fcryst2[i],ffost2[i]]
        chisq_old = chisq[i]
        break
      ENDIF
    ENDFOR

  END
  
  'disk_mips': BEGIN
    fmt = 'a,f,f,f,f,f,f,f,f,f,f,f'
    readcol, 'disk_new.csv',F=fmt,db_name,chisq,rin,rout,rlaw,amin,amax,alaw,diskmass,fcryst,foliv,ffost,/silent

    fmt = 'a,f'
    readcol, 'r_star.csv',F=fmt,r_star_name,c_r_star,/silent

    ; Load correct r_star
    for k=0, n_elements(r_star_name)-1 do begin
      if (r_star_name[k] eq name) then begin
        r_star = c_r_star[k]
        break
      endif
    endfor

    FOR i=0,(n_elements(db_name)-1) DO BEGIN
      IF (name eq db_name[i]) THEN BEGIN
        head = ['rin ','rout ','rlaw ','amin ','amax ','alaw ','diskmass ','folive ','fcrys ','ffors ']
        old_result = [alog10(rin[i]/(r_sun*r_star)),alog(rout[i]/rin[i]),rlaw[i],amin[i],alog(amax[i])/amin[i],alaw[i],alog10((diskmass[i])*m_moon),foliv[i],fcryst[i],ffost[i]]
        chisq_old = chisq[i]
        break
      ENDIF
    ENDFOR

  END
ENDCASE


print, head
print, old_result
print, ' '
print, 'Chisq: ', chisq_old
print, ' '

end