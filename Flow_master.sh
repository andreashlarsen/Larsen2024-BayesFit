

for dir in 2 3 7 8 10 12 14
do

mkdir -p folder$dir
cp Flow_template.sh folder$dir/Flow_folder$dir.sh
cd folder$dir

if [ $dir -eq 2 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0 1 2 3 4 5 6 7 8/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0/g' Flow_folder$dir.sh
elif [ $dir -eq 3 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0 1 2 3 4/g' Flow_folder$dir.sh
elif [ $dir -eq 7 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0 1 2 3 4/g' Flow_folder$dir.sh
sed -i -e 's/M0=400/M0=2000/g' Flow_folder$dir.sh
elif [ $dir -eq 8 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0 1 2 3 4/g' Flow_folder$dir.sh
sed -i -e 's/M0=400/M0=300/g' Flow_folder$dir.sh
sed -i -e 's/M1=50/M1=300/g' Flow_folder$dir.sh
elif [ $dir -eq 10 ]
then
sed -i -e 's/LOOP_F2/for f2 in 0.1 0.2 0.5 1.0 2.0 5.0 10.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0/g' Flow_folder$dir.sh
elif [ $dir -eq 12 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 10 20 11 21 12 22 13 23/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0/g' Flow_folder$dir.sh
elif [ $dir -eq 14 ]
then
sed -i -e 's/LOOP_F2/for f2 in 1.0/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_PRIOR/for j in 0 1 2 3 10 11 12 13 20 21 22 23/g' Flow_folder$dir.sh
sed -i -e 's/LOOP_WEIGHT/for k in 0/g' Flow_folder$dir.sh
sed -i -e 's/KLD=0/KLD=1/g' Flow_folder$dir.sh
fi

cd ..
done
