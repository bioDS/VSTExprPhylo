#enter root directory
root_dir="../VSTExprPhylo"

mkdir BEAST

sh ../../scripts/make_xml.sh stc_zero 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh stc_zero_hvg 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh stc_zero_norm_hvg 6 0.16666667 BDStrictOrdinal $root_dir

sh ../../scripts/make_cont_xml.sh stc_cont_hvg stc_hvg BDStrictContinuous 209  $root_dir





