basename_for_xml=$1
basename=$2
beastdir_basename=$3"_"$basename
minor_dimension=$4
root_dir=$5


python_script_filename=$root_dir"/scripts/continuous2xml.py"
template_filename=$root_dir"/xml_templates/Continuous.xml"
keys_filename="continuous/keys_"$basename".txt"
taxonset_filename="continuous/taxonset_"$basename".txt"
counts_filename="continuous/counts_"$basename".txt"
beastdir_path="BEAST/"$beastdir_basename
xml_filename=$beastdir_path"/"$basename_for_xml".xml"

mkdir $beastdir_path

python3 $python_script_filename $template_filename $taxonset_filename $keys_filename $counts_filename $basename_for_xml $minor_dimension $xml_filename
