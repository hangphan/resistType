while IFS='' read -r line || [[ -n "$line" ]]; do
    #qsub  $1 $line
    #qsub $1 -s $line
    sh $1 $line R00000042
done < "$2"