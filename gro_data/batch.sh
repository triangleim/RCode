Folder_A="/root/setintersect/RCode/data"
for file_a in ${Folder_A}/*; do
    vara=${file_a:30:10}
    echo "$file_a"
    varb="metis-pro-"
    if [ "$vara" = "$varb" ];then
        ../baselines/reorder $file_a -order hybrid
        ../baselines/reorder $file_a -order mloggapa
        ../baselines/reorder $file_a -order metis
        ../baselines/reorder $file_a -order slashburn
        ../baselines/reorder $file_a -order bfsr
        ../baselines/reorder $file_a -order dfs
        # ../baselines/reorder $file_a -order gro
    fi
done    