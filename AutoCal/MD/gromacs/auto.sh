#rm data.txt

total=20000
python generate_data_txt.py $total

for ((line=1; line<=20000; line=line+1))
do
	dataline=$(awk "NR==$line" data.txt)
	Mn=$(echo $dataline | awk '{print $1}')
	Fe=$(echo $dataline | awk '{print $2}')
	Co=$(echo $dataline | awk '{print $3}')
	Ni=$(echo $dataline | awk '{print $4}')
	Zn=$(echo $dataline | awk '{print $5}')
	CN=$(echo $dataline | awk '{print $6}')
	
	gmx_mpi insert-molecules -ci Mn.pdb -nmol $Mn -try 100 -box 3 3 3 -o box.gro
	gmx_mpi insert-molecules -f box.gro -ci Fe.pdb -nmol $Fe -try 100 -o box.gro
	gmx_mpi insert-molecules -f box.gro -ci Co.pdb -nmol $Co -try 100 -o box.gro
	gmx_mpi insert-molecules -f box.gro -ci Ni.pdb -nmol $Ni -try 100 -o box.gro
	gmx_mpi insert-molecules -f box.gro -ci Zn.pdb -nmol $Zn -try 100 -o box.gro
	gmx_mpi insert-molecules -f box.gro -ci CN.pdb -nmol $CN -try 100 -o box.gro

	gmx_mpi editconf -f box.gro -o box.pdb

	rm *.gro.*
	rm *.pdb.*

	python pdb2cif.py box.pdb box.cif
	mkdir $line
	mv box.cif $line
done
