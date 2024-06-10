#enter root directory
root_dir="../VSTExprPhylo"

mkdir BEAST

# selected

sh ../../scripts/make_xml.sh selected_zero 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh selected_zero_hvg 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh selected_zero_norm_hvg 6 0.16666667 BDStrictOrdinal $root_dir

sh ../../scripts/make_cont_xml.sh selected_cont_hvg selected_hvg BDStrictContinous 114 $root_dir


# highest



sh ../../scripts/make_xml.sh t5q_zero 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh t5q_zero_hvg 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh t5q_zero_norm_hvg 6 0.16666667 BDStrictOrdinal $root_dir

sh ../../scripts/make_cont_xml.sh t5q_cont_hvg t5q_hvg BDStrictContinous 114 $root_dir

# similar

sh ../../scripts/make_xml.sh stc_zero 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh stc_zero_hvg 4 0.25 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh stc_zero_norm_hvg 6 0.16666667 BDStrictOrdinal $root_dir

sh ../../scripts/make_cont_xml.sh stc_cont_hvg stc_hvg BDStrictContinous 113  $root_dir





