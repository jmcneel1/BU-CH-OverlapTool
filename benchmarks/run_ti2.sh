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
22
0.0
0.0
0.0
1
22
${nx}
${ny}
${nz}
1
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
2
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
1
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
2
END
echo "YAEHMOP" >> ti2.txt
echo " " >> ti2.txt
cp ../../yaehmop/benchmark_check/ti2_0_0.inp yaehmop/ti2_${i}_${j}.inp
sed -i "s/2\ Ti\ 0.0\ 0.0\ 1.95/2\ Ti\ ${nx}\ ${ny}\ ${nz}/g" yaehmop/ti2_${i}_${j}.inp
cd yaehmop
../../../yaehmop/bind ti2_${i}_${j}.inp
s4s4s="  Ti-4s"
px=" Ti-4px"
py=" Ti-4py"
pz=" Ti-4pz"
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
xypz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n85 | tail -n1 | awk '{print $4}'`
xyx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n85 | tail -n1 | awk '{print $5}'`
xyz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n85 | tail -n1 | awk '{print $6}'`
xzpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n86 | tail -n1 | awk '{print $4}'`
xzx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n86 | tail -n1 | awk '{print $5}'`
xzz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n86 | tail -n1 | awk '{print $6}'`
yzpz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n87 | tail -n1 | awk '{print $4}'`
yzx2y2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n87 | tail -n1 | awk '{print $5}'`
yzz2=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n87 | tail -n1 | awk '{print $6}'`
sxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n98 | tail -n1 | awk '{print $4}'`
sxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n98 | tail -n1 | awk '{print $5}'`
syz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n98 | tail -n1 | awk '{print $6}'`
pxxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n99 | tail -n1 | awk '{print $4}'`
pxxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n99 | tail -n1 | awk '{print $5}'`
pxyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n99 | tail -n1 | awk '{print $6}'`
pyxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n100 | tail -n1 | awk '{print $4}'`
pyxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n100 | tail -n1 | awk '{print $5}'`
pyyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n100 | tail -n1 | awk '{print $6}'`
pzxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n101 | tail -n1 | awk '{print $4}'`
pzxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n101 | tail -n1 | awk '{print $5}'`
pzyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n101 | tail -n1 | awk '{print $6}'`
x2y2xy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n102 | tail -n1 | awk '{print $4}'`
x2y2xz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n102 | tail -n1 | awk '{print $5}'`
x2y2yz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n102 | tail -n1 | awk '{print $6}'`
z2xy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n103 | tail -n1 | awk '{print $4}'`
z2xz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n103 | tail -n1 | awk '{print $5}'`
z2yz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n103 | tail -n1 | awk '{print $6}'`
xyxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n104 | tail -n1 | awk '{print $4}'`
xyxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n104 | tail -n1 | awk '{print $5}'`
xyyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n104 | tail -n1 | awk '{print $6}'`
xzxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n105 | tail -n1 | awk '{print $4}'`
xzxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n105 | tail -n1 | awk '{print $5}'`
xzyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n105 | tail -n1 | awk '{print $6}'`
yzxy=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n106 | tail -n1 | awk '{print $4}'`
yzxz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n106 | tail -n1 | awk '{print $5}'`
yzyz=`grep -A 114 Overlap ti2_${i}_${j}.inp.out | head -n106 | tail -n1 | awk '{print $6}'`
cd ../
echo "         Ti-4s" >> ti2.txt
printf "%s" "${s4s4s}" >> ti2.txt
printf "%7.3f\n\n" "${foursfours}" >> ti2.txt
echo "        Ti-4px Ti-4py Ti-4pz" >> ti2.txt
printf "%s" "${s4s4s}" >> ti2.txt
printf "%7.3f" "${spx}" >> ti2.txt
printf "%7.3f" "${spy}" >> ti2.txt
printf "%7.3f\n\n" "${spz}" >> ti2.txt
echo "              Ti-3xy    Ti-3yz    Ti-3z2    Ti-3xz Ti-3x2y2" >> ti2.txt
printf "     Ti-4s" >> ti2.txt
printf "%10.3f" "${sxy}" >> ti2.txt
printf "%10.3f" "${syz}" >> ti2.txt
printf "%10.3f" "${sz2}" >> ti2.txt
printf "%10.3f" "${sxz}" >> ti2.txt
printf "%10.3f\n\n" "${sx2y2}" >> ti2.txt
echo "         Ti-4s" >> ti2.txt
printf " Ti-4px" >> ti2.txt
printf "%7.3f\n" "${pxs}" >> ti2.txt
printf " Ti-4py" >> ti2.txt
printf "%7.3f\n" "${pys}" >> ti2.txt
printf " Ti-4pz" >> ti2.txt
printf "%7.3f\n\n" "${pzs}" >> ti2.txt
echo "            Ti-4s" >> ti2.txt
printf "   Ti-3dxy" >> ti2.txt
printf "%7.3f\n" "${xys}" >> ti2.txt
printf "   Ti-3dyz" >> ti2.txt
printf "%7.3f\n" "${yzs}" >> ti2.txt
printf "   Ti-3dz2" >> ti2.txt
printf "%7.3f\n" "${z2s}" >> ti2.txt
printf "   Ti-3dxz" >> ti2.txt
printf "%7.3f\n" "${xzs}" >> ti2.txt
printf " Ti-3dx2y2" >> ti2.txt
printf "%7.3f\n\n" "${x2y2s}" >> ti2.txt
echo "        Ti-4px Ti-4py Ti-4pz" >> ti2.txt
printf " Ti-4px" >> ti2.txt
printf "%7.3f" "${pxpx}" >> ti2.txt
printf "%7.3f" "${pxpy}" >> ti2.txt
printf "%7.3f\n" "${pxpz}" >> ti2.txt
printf " Ti-4py" >> ti2.txt
printf "%7.3f" "${pypx}" >> ti2.txt
printf "%7.3f" "${pypy}" >> ti2.txt
printf "%7.3f\n" "${pypz}" >> ti2.txt
printf " Ti-4pz" >> ti2.txt
printf "%7.3f" "${pzpx}" >> ti2.txt
printf "%7.3f" "${pzpy}" >> ti2.txt
printf "%7.3f\n\n" "${pzpz}" >> ti2.txt
echo "          Ti-3dxy   Ti-3dyz   Ti-3dz2   Ti-3dxz Ti-3dx2y2" >> ti2.txt
printf " Ti-4px" >> ti2.txt
printf "%10.3f" "${pxxy}" >> ti2.txt
printf "%10.3f" "${pxyz}" >> ti2.txt
printf "%10.3f" "${pxz2}" >> ti2.txt
printf "%10.3f" "${pxxz}" >> ti2.txt
printf "%10.3f\n" "${pxx2y2}" >> ti2.txt
printf " Ti-4py" >> ti2.txt
printf "%10.3f" "${pyxy}" >> ti2.txt
printf "%10.3f" "${pyyz}" >> ti2.txt
printf "%10.3f" "${pyz2}" >> ti2.txt
printf "%10.3f" "${pyxz}" >> ti2.txt
printf "%10.3f\n" "${pyx2y2}" >> ti2.txt
printf " Ti-4pz" >> ti2.txt
printf "%10.3f" "${pzxy}" >> ti2.txt
printf "%10.3f" "${pzyz}" >> ti2.txt
printf "%10.3f" "${pzz2}" >> ti2.txt
printf "%10.3f" "${pzxz}" >> ti2.txt
printf "%10.3f\n\n" "${pzx2y2}" >> ti2.txt
echo "           Ti-4px Ti-4py Ti-4pz" >> ti2.txt
printf "   Ti-3dxy" >> ti2.txt
printf "%7.3f" "${xypx}" >> ti2.txt
printf "%7.3f" "${xypy}" >> ti2.txt
printf "%7.3f\n" "${xypz}" >> ti2.txt
printf "   Ti-3dyz" >> ti2.txt
printf "%7.3f" "${yzpx}" >> ti2.txt
printf "%7.3f" "${yzpy}" >> ti2.txt
printf "%7.3f\n" "${yzpz}" >> ti2.txt
printf "   Ti-3dz2" >> ti2.txt
printf "%7.3f" "${z2px}" >> ti2.txt
printf "%7.3f" "${z2py}" >> ti2.txt
printf "%7.3f\n" "${z2pz}" >> ti2.txt
printf "   Ti-3dxz" >> ti2.txt
printf "%7.3f" "${xzpx}" >> ti2.txt
printf "%7.3f" "${xzpy}" >> ti2.txt
printf "%7.3f\n" "${xzpz}" >> ti2.txt
printf " Ti-3dx2y2" >> ti2.txt
printf "%7.3f" "${x2y2px}" >> ti2.txt
printf "%7.3f" "${x2y2py}" >> ti2.txt
printf "%7.3f\n\n" "${x2y2pz}" >> ti2.txt
echo "             Ti-3dxy   Ti-3dyz   Ti-3dz2   Ti-3dxz Ti-3dx2y2" >> ti2.txt
printf "   Ti-3dxy" >> ti2.txt
printf "%10.3f" "${xyxy}" >> ti2.txt
printf "%10.3f" "${xyyz}" >> ti2.txt
printf "%10.3f" "${xyz2}" >> ti2.txt
printf "%10.3f" "${xyxz}" >> ti2.txt
printf "%10.3f\n" "${xyx2y2}" >> ti2.txt
printf "   Ti-3dyz" >> ti2.txt
printf "%10.3f" "${yzxy}" >> ti2.txt
printf "%10.3f" "${yzyz}" >> ti2.txt
printf "%10.3f" "${yzz2}" >> ti2.txt
printf "%10.3f" "${yzxz}" >> ti2.txt
printf "%10.3f\n" "${yzx2y2}" >> ti2.txt
printf "   Ti-3dz2" >> ti2.txt
printf "%10.3f" "${z2xy}" >> ti2.txt
printf "%10.3f" "${z2yz}" >> ti2.txt
printf "%10.3f" "${z2z2}" >> ti2.txt
printf "%10.3f" "${z2xz}" >> ti2.txt
printf "%10.3f\n" "${z2x2y2}" >> ti2.txt
printf "   Ti-3dxz" >> ti2.txt
printf "%10.3f" "${xzxy}" >> ti2.txt
printf "%10.3f" "${xzyz}" >> ti2.txt
printf "%10.3f" "${xzz2}" >> ti2.txt
printf "%10.3f" "${xzxz}" >> ti2.txt
printf "%10.3f\n" "${xzx2y2}" >> ti2.txt
printf " Ti-3dx2y2" >> ti2.txt
printf "%10.3f" "${x2y2xy}" >> ti2.txt
printf "%10.3f" "${x2y2yz}" >> ti2.txt
printf "%10.3f" "${x2y2z2}" >> ti2.txt
printf "%10.3f" "${x2y2xz}" >> ti2.txt
printf "%10.3f\n\n" "${x2y2x2y2}" >> ti2.txt
echo " " >> ti2.txt
echo "----------------------------------------------------" >> ti2.txt
done
done
