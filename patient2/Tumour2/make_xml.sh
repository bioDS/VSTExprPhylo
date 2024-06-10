#enter root directory TODO:remove the path
root_dir="/Users/aga122/Git/VSTExprPhylo"

mkdir BEAST

# selected

sh ../../scripts/make_xml.sh testset_zero 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh testset_zero_hvg 5 0.2 BDStrictOrdinal $root_dir

sh ../../scripts/make_xml.sh testset_zero_norm_hvg 6 0.16666667 BDStrictOrdinal $root_dir

sh ../../scripts/make_cont_xml.sh testset_cont_hvg testset_hvg BDStrictContinous 188 $root_dir
