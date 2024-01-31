; .r '/home/jl/jl-data/0_pro_original/sav2vtk/sav2vtk_scalar_qfactor.pro' 

files=file_search('*lev1.fits', count=n)
read_sdo,files[0],header,data
index2map,header,data,aiamap

; QFACTOR
; ===================================================
data=aiamap.data

ss=size(data)
R_SUN=header.R_SUN
P_SIZE=header.cdelt1
c1=header.crpix1
c2=header.crpix2
dy=1/R_SUN
dz=dy
c1_asc=-float(c1)*dy
c2_asc=-float(c2)*dz

;dx=0.25 & dy=0.25 & dz=0.25
;dx=1&dy=1
;xs=ss[1] & ys=ss[2] & zs=ss[3]
filename='aia_211.vtk'
openw,1,filename
printf,1,'# vtk DataFile Version 2.0'
printf,1,'Volume example'
;printf,1,'ASCII'
printf,1,'BINARY'
printf,1,'DATASET STRUCTURED_POINTS'
sn1=strcompress(string(ss[1],format='(i4)'))
sn2=strcompress(string(ss[2],format='(i4)'))
sn=strcompress(string(double(ss[1]*ss[2]),format='(i12)'))
printf,1,'DIMENSIONS '+'1 '+sn1+' '+sn2
;printf,1,'POINTS '+sn+' '+'float'
printf,1,'ORIGIN',0,c1_asc,c2_asc
printf,1,'SPACING',dy,dy,dz
printf,1,'POINT_DATA ',sn
printf,1,'SCALARS decay_index FLOAT'
printf,1,'LOOKUP_TABLE default'
    aia_211=data
    writeu,1,SWAP_ENDIAN(float(aia_211),/swap_if_little_endian)
print,'plane aia is completed !!!'
close,1
end


