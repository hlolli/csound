<CsoundSynthesizer>
<CsInstruments>

sr	=	44100
ksmps	=	1
nchnls	=	2
0dbfs	=	1

opcode testop, 0, 0  
  prints "HELLO WORLD \n"
endop

opcode myadd, i,i 
ival xin

xout ival + 1
endop

opcode testop_new():()
  prints "HELLO WORLD 2\n"
endop

opcode myadd_new(ival):(i)
  prints "In myadd_new\n"
  xout ival + 1
endop


opcode myadd_new2(value0:i):(i)
  xout value0 + 1
endop

instr 1	
 
testop()
testop_new()
ival = myadd(4)
print ival 

ival2 = myadd_new(4)
print ival2 


ival3 = myadd_new2(4)
print ival3 

turnoff

endin

</CsInstruments>
<CsScore>
i1 0 0.2 
</CsScore>
</CsoundSynthesizer>
