stop_on_error ()
## Function that kills the script if you provide it with an input of anything
## except 0 as the first argument. Can be used to interpret error codes,
## since "0" is considered to be a success.
{
    if [ $1 -ne 0 ];
    then
        echo "Last commmand exited with nonzero error code: $1"
        exit
    fi
    return 0
}


logger()
## Logging function. Match arguments and print time.
{
    printf "[@$(date +'%m/%d/%Y %H:%M:%S')]: "
    for i in $*
    do 
        printf "$i " 
    done
    printf "\n"
    return 0
}

q2export()
## Export the file $2 from the qza object $1
{
    term=$2
    inf=$1
    ext=$(echo $term | sed -e 's/.*\.//g')
    outf=$(basename -s .qza $inf).$ext
    zipath=$(unzip -l $inf | grep "$term" | sed 's/.* //g')
    unzip -p $inf $zipath > $outf
}