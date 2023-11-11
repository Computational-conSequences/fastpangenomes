#### mash:
LASTNAME="$1"
if [[ -d fna-${LASTNAME} ]]; then
else
    echo "I need a group name [Enterobacteriaceae]"
    exit 0
fi

mkdir -p Sketch-${LASTNAME} Mash
for FILE in fna-${LASTNAME}/*fna.gz;
do
    SKETCHFL=Sketch-${LASTNAME}/$FILE:t:r:r
    if [[ -f $SKETCHFL.msh ]]; then
    else
        mash sketch -s 5000 -I $FILE:t:r:r -o $SKETCHFL $FILE
    fi
done
#### running mash and have all the results in a table:
MASHTBL="Mash/mash-${LASTNAME}.tbl"
if [ -f $MASHTBL ]; then
    rm $MASHTBL
fi
echo "mashing into $MASHTBL"
mash triangle -p 4 -E Sketch-${LASTNAME}/*.msh > $MASHTBL
bzip2 --best -f $MASHTBL
