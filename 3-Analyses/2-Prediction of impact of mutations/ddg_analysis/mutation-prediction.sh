#!/bin/bash
# Title: mutation-prediction.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-02-29
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

# Compute ddG from wild-type and mutant protein to test impact of point mutation on protein stability



#==========#
# Versions #
#==========#

# v0.0 - 2019-02-29: creation



#===========#
# Functions #
#===========#

# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    if [[ $? != 0 ]]
    then
        error "Package $1 is needed. Exiting..." 1
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 1
}



#==============#
# Dependencies #
#==============#

test_dep obabel
test_dep minimize_with_cst.static.linuxgccrelease
test_dep ddg_monomer.static.linuxgccrelease



#===========#
# Variables #
#===========#

myfilename="Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel"

rosetta_path="$HOME/local/pckg/rosetta_bin_linux_3.9_bundle/main"
[[ ! -d "$rosetta_path" ]] && error "Cannot detect Rosetta path. Exiting..." 1



#============#
# Processing #
#============#

#-------------#
# Param files #
#-------------#

info "Generating param files..."

cd params

# Add hydrogens
obabel OAQ.sdf -O OAQ_withH.sdf -h
obabel A3P.sdf -O A3P_withH.sdf -h

# Generating params files
"$rosetta_path/source/scripts/python/public/molfile_to_params.py" -n OAQ -p OAQ --conformers-in-one-file OAQ_withH.sdf
"$rosetta_path/source/scripts/python/public/molfile_to_params.py" -n A3P -p A3P --conformers-in-one-file A3P_withH.sdf

cd ..


#-----------------#
# Constraint file #
#-----------------#

info "Generating constraint file..."

# Create a list (-in:file:s gives error)
echo 4MUB.pdb > lst

# Minimize structure
minimize_with_cst.static.linuxgccrelease -in:file:l lst \
    -in:file:fullatom \
    -extra_res_fa params/OAQ.params \
    -extra_res_fa params/A3P.params \
    -fa_max_dis 9.0 \
    -database "$rosetta_path/database/" \
    -ddg::harmonic_ca_tether 0.5     \
    -score:weights score12           \
    -ddg::constraint_weight 1.0      \
    -ddg::out_pdb_prefix min_cst_0.5 \
    -ddg::sc_min_only false          \
    -score:patch score12.wts_patch > mincst.log

# Convert log in cst (if run for the first time, script might required a chmod u+x)
"$rosetta_path/source/src/apps/public/ddg/convert_to_cst_file.sh" mincst.log > 4MUB.cst

# Update the cst file to remove A3P and OAQ information
egrep -wv "257|258" 4MUB.cst > 4MUB.cst.tmp
mv 4MUB.cst.tmp 4MUB.cst


#-----------------#
# ddG computation #
#-----------------#

info "Computing ddG (this will take a while)..."

# List of mutations. Run each line independently if run on several jobs / servers in parallel
mutations="A74T C35R D152E \
    D174N D223A E142K \
    E180K E3G G204D \
    G204S G206V H122N \
    H251Y H37N H57N \
    I50M K244Q K71T \
    K99N L179P L188I \
    L188M L256W M252I \
    N215Y N44D P106S \
    P130H P209S P225S \
    Q176R R202H S160L \
    T89A V189A V189I \
    V227I W120R W52C Y243H"

# Iteration (can be whatever number but the higher it is, the longer it runs)
iteration=50
myit_dir="iteration_${iteration}_pre-talaris"

for i in $mutations
do
    echo "$i"

    # Create folder or skip if already existing
    [[ -d "$myit_dir/$i" ]] && echo "Folder exists, skipping..." && continue || mkdir -p "$myit_dir/$i"
    
    cd "$myit_dir/$i"
    
    # Link needed files
    ln -s ../../min_cst_0.5.4MUB_0001.pdb .
    ln -s ../../4MUB.cst .
    ln -s ../../params .
    
    # Create resfile
    echo "\
NATAA
USE_INPUT_SC    # allow the use of the input side chain conformation
start
$(echo "$i" | sed "s/[A-Z]//g") A PIKAA $(echo "$i" | sed -r "s/.*[0-9].*([A-Z])/\1/g")" > resfile
    
    # Compute ddG with defined number of iterations
    ddg_monomer.static.linuxgccrelease -in:file:s min_cst_0.5.4MUB_0001.pdb \
        -resfile resfile \
        -ddg:weight_file soft_rep_design        \
        -ddg:minimization_scorefunction score12 \
        -database "$rosetta_path/database/"\
        -fa_max_dis 9.0              \
        -ddg::iterations $iteration  \
        -ddg::dump_pdbs true         \
        -ddg::local_opt_only false   \
        -ddg::min_cst true           \
        -constraints::cst_file 4MUB.cst   \
        -ddg::suppress_checkpointing true \
        -in::file::fullatom \
        -ddg::mean false    \
        -ddg::min true      \
        -ddg::sc_min_only false   \
        -ddg::ramp_repulsive true \
        -mute all                 \
        -unmute core.optimization.LineMinimizer \
        -ddg::output_silent true        \
        -extra_res_fa params/OAQ.params \
        -extra_res_fa params/A3P.params \
        -restore_pre_talaris_2013_behavior true &> ddg_monomer.log # Restoring behavior is critical!
    
    cd ../..
done


#-----------#
# ddG table #
#-----------#

info "Generating ddG table..."

for i in $(ls -1 iteration_50_pre-talaris/)
do
    echo -e "$i\t$(cat iteration_50_pre-talaris/$i/ddg_p* | sed "s/  */\t/g" | cut -f 3 | egrep [0-9] | tail -1)"
done > ddg_mutations.tsv

exit 0
