mol load psf step5_charmm2namd.pdb xtc final1_5us.xtc
set nf [molinfo top get numframes]

set sel1 [ atomselect top "protein and resid 214 and name CA"]
set sel2 [ atomselect top "protein and 338 and name CA"]

set outfile [open "Ionic_lock.dat" w]

for {set i 0} {$i < $nf} {incr i} {
puts "frame $i of $nf"
$sel1 frame $i
$sel2 frame $i
set com1 [measure center $sel1 weight mass]
set com2 [measure center $sel2 weight mass]
set simdata($i.r) [veclength [vecsub $com1 $com2]]
puts $outfile "$i , $simdata($i.r)"
}

close $outfile


#set sel1 [ atomselect top "protein and resid 101 and name CE1 CE2 CZ CD1 CD2 CG and noh"]
#set sel2 [ atomselect top "protein and resid 257 and name CG CD1 NE1 CD2 CE2 CZ2 CH2 CZ3 CE3 and noh"]




