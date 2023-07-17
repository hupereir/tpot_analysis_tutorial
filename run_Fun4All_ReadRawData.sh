#!/bin/csh
source /opt/sphenix/core/bin/sphenix_setup.csh 

set nEvents = ${1}
set runNumber = ${2}
set type = ${3}

# input file name
set input_path = "/sphenix/lustre01/sphnxpro/commissioning/TPOT/${type}/"
set inputfile = `printf "${input_path}/TPOT_ebdc39_${type}-%08i-0000.prdf" $runNumber`
echo "inputfile: ${inputfile}"

# check input file existence
if ( ! -f ${inputfile} ) then
  echo "File ${inputfile} does not exists. aborting."
  exit 0
endif

# output path and filename
set output_path = "CONDOR_RawDataEvaluation/"
set outputfile = `printf "${output_path}/MicromegasRawDataEvaluation-%08i-0000.root" $runNumber`
echo "outputfile: ${outputfile}"

# make sure output path exists
mkdir -p ${output_path}

# check output file existence
if ( -f ${outputfile} ) then
  echo "File ${outputfile} exists. aborting."
  exit 0
endif

touch ${outputfile}

# start root, load macro, run and exit
echo "running root"
root -b << EOF
.L Fun4All_ReadRawData.C
Fun4All_ReadRawData( $nEvents, "${inputfile}", "${outputfile}" )
EOF
