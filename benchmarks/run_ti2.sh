for i in `seq 0 30 330`; do 
for j in `seq 0 30 330`; do
echo " " >> ti2.txt
echo "r = 1.95 A, theta = ${i}, phi = ${j}" >> ti2.txt
echo " " >> ti2.txt 
echo "Overlap_Tool" >> ti2.txt
echo " " >> ti2.txt
nx=`echo "1.95*s(${i}/180*3.14159265359)*c(${j}/180*3.14159265359)" | bc -l` 
ny=`echo "1.95*s(${i}/180*3.14159265359)*s(${j}/180*3.14159265359)" | bc -l`
nz=`echo "1.95*c(${i}/180*3.14159265359)" | bc -l`
../overlap_tool << END | tail -n3 >> ti2.txt
22
0.0
0.0
0.0
0
22
${nx}
${ny}
${nz}
0
END
../overlap_tool << END | tail -n3 >> ti2.txt
22
0.0
0.0
0.0
0
22
${nx}
${ny}
${nz}
1
END
../overlap_tool << END | tail -n3 >> ti2.txt
22
0.0
0.0
0.0
0
22
${nx}
${ny}
${nz}
2
END
../overlap_tool << END | tail -n5 >> ti2.txt
22
0.0
0.0
0.0
1
22
${nx}
${ny}
${nz}
0
END
../overlap_tool << END | tail -n7 >> ti2.txt
22
0.0
0.0
0.0
2
22
${nx}
${ny}
${nz}
0
END
../overlap_tool << END | tail -n5 >> ti2.txt
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
echo "YAEHMOP" >> ti2.txt
echo " " >> ti2.txt
cp ../../yaehmop/benchmark/ti2_0_0.inp yaehmop/ti2_${i}_${j}.inp
sed -i "s/2\ Ti\ 0.0\ 0.0\ 1.95/2\ Ti\ ${nx}\ ${ny}\ ${nz}/g" yaehmop/ti2_${i}_${j}.inp
cd yaehmop
../../../yaehmop/bind ti2_${i}_${j}.inp
s4s4s="  Ti-2s"
spx=" Ti-2px"
spy=" Ti-2py"
spz=" Ti-2pz"
foursfours=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n60 | tail -n1 | awk '{print $4}'`
spx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n60 | tail -n1 | awk '{print $5}'`
spy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n60 | tail -n1 | awk '{print $6}'`
pxs=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n61 | tail -n1 | awk '{print $4}'`
pxpx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n61 | tail -n1 | awk '{print $5}'`
pxpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n61 | tail -n1 | awk '{print $6}'`
pys=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n62 | tail -n1 | awk '{print $4}'`
pypx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n62 | tail -n1 | awk '{print $5}'`
pypy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n62 | tail -n1 | awk '{print $6}'`
pzs=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n63 | tail -n1 | awk '{print $4}'`
pzpx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n63 | tail -n1 | awk '{print $5}'`
pzpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n63 | tail -n1 | awk '{print $6}'`
x2y2s=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n64 | tail -n1 | awk '{print $4}'`
x2y2px=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n64 | tail -n1 | awk '{print $5}'`
x2y2py=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n64 | tail -n1 | awk '{print $6}'`
z2s=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n65 | tail -n1 | awk '{print $4}'`
z2px=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n65 | tail -n1 | awk '{print $5}'`
z2py=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n65 | tail -n1 | awk '{print $6}'`
xys=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n66 | tail -n1 | awk '{print $4}'`
xypx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n66 | tail -n1 | awk '{print $5}'`
xypy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n66 | tail -n1 | awk '{print $6}'`
xzs=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n67 | tail -n1 | awk '{print $4}'`
xzpx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n67 | tail -n1 | awk '{print $5}'`
xzpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n67 | tail -n1 | awk '{print $6}'`
yzs=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n68 | tail -n1 | awk '{print $4}'`
yzpx=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n68 | tail -n1 | awk '{print $5}'`
yzpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n68 | tail -n1 | awk '{print $6}'`
spz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $4}'`
sx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $5}'`
sz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $6}'`
pxpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n80 | tail -n1 | awk '{print $4}'`
pxx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n80 | tail -n1 | awk '{print $5}'`
pxz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n80 | tail -n1 | awk '{print $6}'`
pypz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n81 | tail -n1 | awk '{print $4}'`
pyx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n81 | tail -n1 | awk '{print $5}'`
pyz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n81 | tail -n1 | awk '{print $6}'`
pzpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n82 | tail -n1 | awk '{print $4}'`
pzx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n82 | tail -n1 | awk '{print $5}'`
pzz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n82 | tail -n1 | awk '{print $6}'`
x2y2pz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n83 | tail -n1 | awk '{print $4}'`
x2y2x2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n83 | tail -n1 | awk '{print $5}'`
x2y2z2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n83 | tail -n1 | awk '{print $6}'`
z2pz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n84 | tail -n1 | awk '{print $4}'`
z2x2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n84 | tail -n1 | awk '{print $5}'`
z2z2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n84 | tail -n1 | awk '{print $6}'`
spz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $4}'`
sx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $5}'`
sz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $6}'`
spz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $4}'`
sx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $5}'`
sz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $6}'`
spz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $4}'`
sx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $5}'`
sz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n79 | tail -n1 | awk '{print $6}'`
spz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n21 | tail -n1 | awk '{print $5}'`
pxpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n22 | tail -n1 | awk '{print $4}'`
pxpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n22 | tail -n1 | awk '{print $5}'`
pypy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n23 | tail -n1 | awk '{print $4}'`
pypz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n23 | tail -n1 | awk '{print $5}'`
pzpy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n24 | tail -n1 | awk '{print $4}'`
pzpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n24 | tail -n1 | awk '{print $5}'`
cd ../
echo "          F-2s" >> ti2.txt
printf "%s" "${s2s2s}" >> ti2.txt
printf "%7.3f\n\n" "${twostwos}" >> ti2.txt
echo "         F-2px  F-2py  F-2pz" >> ti2.txt
printf "%s" "${s2s2s}" >> ti2.txt
printf "%7.3f" "${spx}" >> ti2.txt
printf "%7.3f" "${spy}" >> ti2.txt
printf "%7.3f\n\n" "${spz}" >> ti2.txt
echo "          F-2s" >> ti2.txt
printf "  F-2px" >> ti2.txt
printf "%7.3f\n" "${pxs}" >> ti2.txt
printf "  F-2py" >> ti2.txt
printf "%7.3f\n" "${pys}" >> ti2.txt
printf "  F-2pz" >> ti2.txt
printf "%7.3f\n\n" "${pzs}" >> ti2.txt
echo "         F-2px  F-2py  F-2pz" >> ti2.txt
printf "  F-2px" >> ti2.txt
printf "%7.3f" "${pxpx}" >> ti2.txt
printf "%7.3f" "${pxpy}" >> ti2.txt
printf "%7.3f\n" "${pxpz}" >> ti2.txt
printf "  F-2py" >> ti2.txt
printf "%7.3f" "${pypx}" >> ti2.txt
printf "%7.3f" "${pypy}" >> ti2.txt
printf "%7.3f\n" "${pypz}" >> ti2.txt
printf "  F-2pz" >> ti2.txt
printf "%7.3f" "${pzpx}" >> ti2.txt
printf "%7.3f" "${pzpy}" >> ti2.txt
printf "%7.3f\n" "${pzpz}" >> ti2.txt
echo " " >> ti2.txt
echo "----------------------------------------------------" >> ti2.txt
done
done
