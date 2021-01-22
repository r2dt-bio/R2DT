for file in *.ps
do
 pstopdf "$file" || echo "problem with $file"
done
