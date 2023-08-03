#!/bin/bash
input="cif_files/cif_library.txt"
while IFS= read -r line
do
    #   echo "$line"
    mkdir $line
    cd $line
    copy every file with the name $line to the new directory
    cp ../cif_files/$line*.cif .
    
    # Make new file and write to it
    cat > $line.py << EOF

    ENTER PYTHON CODE HERE

    EOF
    
    # Run python script but first go into screen session

    cd ..
    
done < "$input"