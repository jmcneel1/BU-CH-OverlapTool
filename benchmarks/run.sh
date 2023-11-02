for i in `seq 0 30 330`; do 
for j in `seq 0 30 330`; do
echo " " >> f2_new.txt
echo "r = 1.45 A, theta = ${i}, phi = ${j}" >> f2_new.txt
echo " " >> f2_new.txt 
echo "Overlap_Tool" >> f2_new.txt
echo " " >> f2_new.txt
nx=`echo "1.45*s(${i}/180*3.14159265359)*c(${j}/180*3.14159265359)" | bc -l` 
ny=`echo "1.45*s(${i}/180*3.14159265359)*s(${j}/180*3.14159265359)" | bc -l`
nz=`echo "1.45*c(${i}/180*3.14159265359)" | bc -l`
../overlap_tool << END | tail -n3 >> f2_new.txt
9
0.0
0.0
0.0
0
9
${nx}
${ny}
${nz}
0
END
../overlap_tool << END | tail -n3 >> f2_new.txt
9
0.0
0.0
0.0
0
9
${nx}
${ny}
${nz}
1
END
../overlap_tool << END | tail -n5 >> f2_new.txt
9
0.0
0.0
0.0
1
9
${nx}
${ny}
${nz}
0
END
../overlap_tool << END | tail -n5 >> f2_new.txt
9
0.0
0.0
0.0
1
9
${nx}
${ny}
${nz}
1
END
echo "YAEHMOP" >> f2_new.txt
echo " " >> f2_new.txt
cp ../../yaehmop/benchmark/f2_0_0.inp yaehmop/f2_${i}_${j}.inp
sed -i "s/2\ F\ 0.0\ 0.0\ 1.45/2\ F\ ${nx}\ ${ny}\ ${nz}/g" yaehmop/f2_${i}_${j}.inp
cd yaehmop
../../../yaehmop/bind f2_${i}_${j}.inp
s2s2s="   F-2s"
spx="  F-2px"
spy="  F-2py"
spz="  F-2pz"
twostwos=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n12 | tail -n1 | awk '{print $5}'`
spx=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n12 | tail -n1 | awk '{print $6}'`
pxs=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n13 | tail -n1 | awk '{print $5}'`
pxpx=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n13 | tail -n1 | awk '{print $6}'`
pys=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n14 | tail -n1 | awk '{print $5}'`
pypx=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n14 | tail -n1 | awk '{print $6}'`
pzs=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n15 | tail -n1 | awk '{print $5}'`
pzpx=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n15 | tail -n1 | awk '{print $6}'`
spy=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n21 | tail -n1 | awk '{print $4}'`
spz=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n21 | tail -n1 | awk '{print $5}'`
pxpy=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n22 | tail -n1 | awk '{print $4}'`
pxpz=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n22 | tail -n1 | awk '{print $5}'`
pypy=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n23 | tail -n1 | awk '{print $4}'`
pypz=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n23 | tail -n1 | awk '{print $5}'`
pzpy=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n24 | tail -n1 | awk '{print $4}'`
pzpz=`grep -A 27 Overlap f2_${i}_${j}.inp.out | head -n24 | tail -n1 | awk '{print $5}'`
cd ../
echo "          F-2s" >> f2_new.txt
printf "%s" "${s2s2s}" >> f2_new.txt
printf "%7.3f\n\n" "${twostwos}" >> f2_new.txt
echo "         F-2px  F-2py  F-2pz" >> f2_new.txt
printf "%s" "${s2s2s}" >> f2_new.txt
printf "%7.3f" "${spx}" >> f2_new.txt
printf "%7.3f" "${spy}" >> f2_new.txt
printf "%7.3f\n\n" "${spz}" >> f2_new.txt
echo "          F-2s" >> f2_new.txt
printf "  F-2px" >> f2_new.txt
printf "%7.3f\n" "${pxs}" >> f2_new.txt
printf "  F-2py" >> f2_new.txt
printf "%7.3f\n" "${pys}" >> f2_new.txt
printf "  F-2pz" >> f2_new.txt
printf "%7.3f\n\n" "${pzs}" >> f2_new.txt
echo "         F-2px  F-2py  F-2pz" >> f2_new.txt
printf "  F-2px" >> f2_new.txt
printf "%7.3f" "${pxpx}" >> f2_new.txt
printf "%7.3f" "${pxpy}" >> f2_new.txt
printf "%7.3f\n" "${pxpz}" >> f2_new.txt
printf "  F-2py" >> f2_new.txt
printf "%7.3f" "${pypx}" >> f2_new.txt
printf "%7.3f" "${pypy}" >> f2_new.txt
printf "%7.3f\n" "${pypz}" >> f2_new.txt
printf "  F-2pz" >> f2_new.txt
printf "%7.3f" "${pzpx}" >> f2_new.txt
printf "%7.3f" "${pzpy}" >> f2_new.txt
printf "%7.3f\n" "${pzpz}" >> f2_new.txt
echo " " >> f2_new.txt
echo "----------------------------------------------------" >> f2_new.txt
done
done
