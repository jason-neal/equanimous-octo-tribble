# Bulk rename the badly labled crires files.
for i in $(ls CRIRES*:*:*.fits); do
    python ~/bin/file_renamer.py $i -x : -v -
done
