nphi=360L
ntheta=180L
nrho=192L

bx=fltarr(nphi*(ntheta+1)*(nrho+1))
path='bx_synoptic.dat'
bx_true=fltarr(nphi,ntheta+1,nrho+1)
openr, lun, path, /get_lun
readf,lun,bx
close,lun
free_lun,lun
for i=0,nphi-1 do begin
for j=0,ntheta do begin
for k=0,nrho do begin
bx_true[i,j,k]=bx[i*(ntheta+1)*(nrho+1)+j*(nrho+1)+k]
endfor
endfor
endfor

by=fltarr(nphi*(ntheta+1)*(nrho+1))
path='by_synoptic.dat'
by_true=fltarr(nphi,ntheta+1,nrho+1)
openr, lun, path, /get_lun
readf,lun,by
close,lun
free_lun,lun
for i=0,nphi-1 do begin
for j=0,ntheta do begin
for k=0,nrho do begin
by_true[i,j,k]=by[i*(ntheta+1)*(nrho+1)+j*(nrho+1)+k]
endfor
endfor
endfor

bz=fltarr(nphi*(ntheta+1)*(nrho+1))
path='bz_synoptic.dat'
bz_true=fltarr(nphi,ntheta+1,nrho+1)
openr, lun, path, /get_lun
readf,lun,bz
close,lun
free_lun,lun
for i=0,nphi-1 do begin
for j=0,ntheta do begin
for k=0,nrho do begin
bz_true[i,j,k]=bz[i*(ntheta+1)*(nrho+1)+j*(nrho+1)+k]
endfor
endfor
endfor

rhos=findgen(nrho+1)
rhos=(rhos/nrho*(2.5-1.))+1.
lat=asin((findgen(ntheta+1)-90.)*0.011111111111111112)/!dpi*180.
lon=findgen(nphi+1)-61.6
solar_B=-7.071*!dpi/180.
pi=!dpi
x1c=rhos;161
x2c=0.5*!dpi-((lat)*!dpi/180);385 seta
x3c=lon*!dpi/180;513 fine

b1c=bx_true
b2c=by_true
b3c=bz_true

dim=size(bz_true)
dimx = dim(1);513
dimy = dim(2);385
dimz = dim(3);161

XX=fltarr(dimx,dimy,dimz)
YY=fltarr(dimx,dimy,dimz)
ZZ=fltarr(dimx,dimy,dimz)
BX=fltarr(dimx,dimy,dimz)
BY=fltarr(dimx,dimy,dimz)
BZ=fltarr(dimx,dimy,dimz)

x_i=0.0
y_i=0.0
z_i=0.0
Bx_i=0.0
By_i=0.0
Bz_i=0.0
for k = 0, dimz-1 do begin
  print,'k= ', k
  for j = 0, dimy-1 do begin
    for i = 0, dimx-1 do begin

      x_i=(x1c(k)*sin(x2c(j))*cos(x3c(i)))
      y_i=x1c(k)*sin(x2c(j))*sin(x3c(i))
      z_i=x1c(k)*cos(x2c(j))

      XX(i,j,k)=x_i*cos(-solar_B)-z_i*sin(-solar_B)
      YY(i,j,k)=y_i
      ZZ(i,j,k)=x_i*sin(-solar_B)+z_i*cos(-solar_B)
      
      bx_i=(-b1c(i,j,k)*sin(x3c(i))-$
                b2c(i,j,k)*cos(x2c(j))*cos(x3c(i))+$
                b3c(i,j,k)*sin(x2c(j))*cos(x3c(i)))
      by_i=b1c(i,j,k)*cos(x3c(i))-$
                b2c(i,j,k)*cos(x2c(j))*sin(x3c(i))+$
                b3c(i,j,k)*sin(x2c(j))*sin(x3c(i))
      bz_i=(b2c(i,j,k)*sin(x2c(j))+b3c(i,j,k)*cos(x2c(j)))

      BX(i,j,k)=bx_i*cos(-solar_B)-bz_i*sin(-solar_B)
      BY(i,j,k)=by_i
      BZ(i,j,k)=bx_i*sin(-solar_B)+bz_i*cos(-solar_B)
      
    endfor
  endfor
endfor
arry=b1c


ss=size(arry)
print,ss
;dx=0.25 & dy=0.25 & dz=0.25
;xs=ss[1] & ys=ss[2] & zs=ss[3]
filename='Synoptic.vtk'
openw,1,filename
printf,1,'# vtk DataFile Version 2.0'
printf,1,'Volume example'
;printf,1,'ASCII'
printf,1,'BINARY'
printf,1,'DATASET STRUCTURED_GRID'
sn1=strcompress(string(ss[1],format='(i4)'))
sn2=strcompress(string(ss[2],format='(i4)'))
sn3=strcompress(string(ss[3],format='(i4)'))
sn=strcompress(string(double(ss[1]*ss[2]*ss[3]),format='(i12)'))
printf,1,'DIMENSIONS '+sn1+' '+sn2+' '+sn3
printf,1,'POINTS '+sn+' '+'float'
for k=0,ss[3]-1 do begin
    for j=0,ss[2]-1 do begin
        for i=0,ss[1]-1 do begin
        ;printf,1,xx[i,j,k],yy[i,j,k],zz[i,j,k],format='(3E18.8)'
	vector_l=[xx[i,j,k],yy[i,j,k],zz[i,j,k]]
	writeu,1,SWAP_ENDIAN(float(vector_l),/swap_if_little_endian)
        endfor
   endfor
endfor

print,'Lvector field is completed !!!'
printf,1,'POINT_DATA', sn
;printf,1,'ORIGIN',0,0,0
;printf,1,'SPACING',dx,dy,dz
;printf,1,'POINT_DATA ',sn
;printf,1,'SCALARS Scalars_twist FLOAT'
;printf,1,'LOOKUP_TABLE default'
printf,1,'VECTORS Vector_B FLOAT'
;printf,1,'LOOKUP_TABLE default'
for k=0,ss[3]-1 do begin
    for j=0,ss[2]-1 do begin
      for i=0,ss[1]-1 do begin
        vector_b=[bx[i,j,k],by[i,j,k],bz[i,j,k]]
	writeu,1,SWAP_ENDIAN(float(vector_b),/swap_if_little_endian)
       endfor
  endfor
endfor    
print,'Bvector field is completed !!!'
close,1
end
