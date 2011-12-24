'reinit'
'open date.ctl'
'enable print love.gmf'
i=1
while(i<25)
'draw title the 'i 'hour U,V,Z filed'
'set t 'i
'set gxout stream'
'set cint 40'
'd z'
'd u;v'
'print'
i=i+1
'c'
endwhile
'disable print'
'enable print Height.gmf'
i=1
while(i<25)
'draw title the 'i 'hour Height filed'
'set t 'i
'set cint 20'
'd z'
'print'
i=i+1
'c'
endwhile
'disable print'
'reinit'