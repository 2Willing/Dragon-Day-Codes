;cd,'/home/test/modelsav/'
file = './model7_30000_export.sav'

; TWIST
; ===================================================

restore,file
;b1c=b[*,*,*,0]
;b2c=b[*,*,*,1]
b3c=b[*,*,*,2]

;b1c=reform(b1,161,385,513)
;b2c=reform(b2,161,385,513)
;b3c=reform(b3,161,385,513)

solar_B=-7.071*!dpi/180.
pi=!dpi
x1c=rad;161
x2c=0.5*!dpi-((lat)*!dpi/180);385 seta
x3c=lon*!dpi/180;513 fine

dim = size(b3c)
dimx = dim(1);513
dimy = dim(2);385
dimz = 1;dim(3);161

XX=fltarr(dimx,dimy);,dimz)
YY=fltarr(dimx,dimy);,dimz)
ZZ=fltarr(dimx,dimy);,dimz)
BZ=fltarr(dimx,dimy);,dimz)


x_i=0.0
y_i=0.0
z_i=0.0
k=0
print,'k= ', k
for j = 0, dimy-1 do begin
  for i = 0, dimx-1 do begin

    x_i=(x1c(k)*sin(x2c(j))*cos(x3c(i)))
    y_i=x1c(k)*sin(x2c(j))*sin(x3c(i))
    z_i=x1c(k)*cos(x2c(j))

    XX(i,j)=x_i*cos(-solar_B)-z_i*sin(-solar_B)
    YY(i,j)=y_i
    ZZ(i,j)=x_i*sin(-solar_B)+z_i*cos(-solar_B)
    
    BZ(i,j)=b3c(i,j,k);*sin(x2c(j))+b3c(i,j,k)*cos(x2c(j)))

    ;BX(i,j,k)=bx_i*cos(-solar_B)-bz_i*sin(-solar_B)
    ;BY(i,j,k)=by_i
    ;BZ(i,j,k)=bx_i*sin(-solar_B)+bz_i*cos(-solar_B)
    
  endfor
endfor
arry=bz


ss=size(arry)
print,ss
;dx=0.25 & dy=0.25 & dz=0.25
;xs=ss[1] & ys=ss[2] & zs=ss[3]
filename='Bz_sphere.vtk'
openw,1,filename
printf,1,'# vtk DataFile Version 2.0'
printf,1,'Volume example'
;printf,1,'ASCII'
printf,1,'BINARY'
printf,1,'DATASET STRUCTURED_GRID'
sn1=strcompress(string(ss[1],format='(i4)'))
sn2=strcompress(string(ss[2],format='(i4)'))
sn3=strcompress(string(1,format='(i4)'))
sn=strcompress(string(double(ss[1]*ss[2]),format='(i12)'))
printf,1,'DIMENSIONS '+sn1+' '+sn2+' '+sn3
printf,1,'POINTS '+sn+' '+'float'
k=0
  for j=0,ss[2]-1 do begin
      for i=0,ss[1]-1 do begin
      ;printf,1,xx[i,j,k],yy[i,j,k],zz[i,j,k],format='(3E18.8)'
vector_l=[xx[i,j],yy[i,j],zz[i,j]]
writeu,1,SWAP_ENDIAN(float(vector_l),/swap_if_little_endian)
      endfor
  endfor

print,'Lvector field is completed !!!'
printf,1,'POINT_DATA', sn
;printf,1,'ORIGIN',0,0,0
;printf,1,'SPACING',dx,dy,dz
;printf,1,'POINT_DATA ',sn
;printf,1,'SCALARS Scalars_twist FLOAT'
;printf,1,'LOOKUP_TABLE default'
printf,1,'SCALARS B_z FLOAT'
printf,1,'LOOKUP_TABLE default'
k=0
    for j=0,ss[2]-1 do begin
      for i=0,ss[1]-1 do begin
        bz_i=bz[i,j]
	writeu,1,SWAP_ENDIAN(float(bz_i),/swap_if_little_endian)
       endfor
  endfor 
print,'Bz field is completed !!!'
close,1
end
