# download jar for flagger
wget -O ./bin/VolcanoSV-asm/Flagger/cromwell-85.jar https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
wget -O ./bin/VolcanoSV-asm/Flagger/womtool-85.jar https://github.com/broadinstitute/cromwell/releases/download/85/womtool-85.jar

# build miniasm
cd ./bin/VolcanoSV-asm/miniasm;make
