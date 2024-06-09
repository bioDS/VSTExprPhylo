fasta_basename=$1
stateNr=$2
frequency=$3
analyses=$4
beastdir_basename=$4"_"$1
root_dir=$5

python_script_filename=$root_dir"/scripts/fasta2xml_expr.py"
template_filename=$root_dir"/xml_templates/BDStrictOrdinal_zero.xml"
fasta_filename="fasta/"$fasta_basename".fasta"
beastdir_path="BEAST/"$beastdir_basename
xml_filename=$beastdir_path"/"$fasta_basename".xml"

mkdir $beastdir_path

python3 $python_script_filename $template_filename $fasta_filename $xml_filename

sed 's/numberofstates/'$stateNr'/g' $xml_filename > tmp
sed 's/initialfrequency/'$frequency'/g' tmp > $xml_filename
sed 's/logfilename/'$fasta_basename'/g' $xml_filename > tmp
sed 's/treefilename/'$fasta_basename'/g' tmp > $xml_filename

rm tmp
