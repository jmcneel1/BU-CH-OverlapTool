for i in `seq 0 30 330`; do 
echo "${i}"
nx=`echo "1.45*s(${i}/180*3.14159265359)" | bc -l` 
nz=`echo "1.45*c(${i}/180*3.14159265359)" | bc -l`
./overlap_tool << END
9
0.0
0.0
0.0
1
9
${nx}
0.0
${nz}
1

END
done
